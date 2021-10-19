#
# DCore -- Integrated DMFT software for correlated electrons
# Copyright (C) 2017 The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#


import os
import sys
import copy
import numpy
import h5py
import ast
import time

from h5 import HDFArchive

from dcore.dmft_core import DMFTCoreSolver
from dcore.program_options import create_parser, parse_parameters
from dcore.tools import *
from dcore import impurity_solvers
from .sumkdft_workers.launcher import run_sumkdft

# from BSE repo
from bse_tools.h5bse import h5BSE


def to_str(x):
    if isinstance(x, bytes):
        return x.decode('utf-8')
    return x


def compare_str_list(list1, list2):
    if len(list1) != len(list2):
        return False
    for x, y in zip(list1, list2):
        if to_str(x) != to_str(y):
            return False
    return True


def calc_g2_in_impurity_model(solver_name, solver_params, mpirun_command, basis_rot, Umat, gf_struct, beta, n_iw,
                              Sigma_iw, Gloc_iw, num_wb, num_wf, ish, freqs=None):
    """

    Calculate G2 in an impurity model

    if freqs is None: box frequency sampling
    else: sparse sampling

    """

    Solver = impurity_solvers.solver_classes[solver_name]

    raise_if_mpi_imported()

    sol = Solver(beta, gf_struct, Umat, n_iw)

    G0_iw = dyson(Sigma_iw=Sigma_iw, G_iw=Gloc_iw)
    sol.set_G0_iw(G0_iw)

    # Compute rotation matrix to the diagonal basis if supported
    rot = compute_diag_basis(G0_iw) if basis_rot else None
    s_params = copy.deepcopy(solver_params)
    s_params['random_seed_offset'] = 1000 * ish

    work_dir_org = os.getcwd()
    work_dir = 'work/imp_shell' + str(ish) + '_bse'
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)

    # True: box sampling,  False: sparse sampling
    flag_box = freqs is None

    # Solve the model
    if flag_box:
        xloc, chiloc = sol.calc_Xloc_ph(rot, mpirun_command, num_wf, num_wb, s_params)
    else:
        xloc, chiloc = sol.calc_Xloc_ph_sparse(rot, mpirun_command, freqs, num_wb, s_params)

    # Check results for x_loc
    print("\n Checking x_loc...")
    assert isinstance(xloc, dict)
    for key, data in list(xloc.items()):
        # print("  ", key)
        if flag_box:
            assert data.shape == (num_wb, 2*num_wf, 2*num_wf)
        else:
            assert data.shape == (freqs.shape[0],)
    print(" OK")

    # Check results for chi_loc
    if chiloc is not None:
        print("\n Checking chi_loc...")
        assert isinstance(chiloc, dict)
        for key, data in list(chiloc.items()):
            # print("  ", key)
            assert data.shape == (num_wb, )
        print(" OK")

    os.chdir(work_dir_org)

    return xloc, chiloc, sol.get_Gimp_iw()


def subtract_disconnected(xloc, gimp, spin_names, freqs=None):
    """
    subtract disconnected part from X_loc
        X_loc[(i1, i2, i3, i4)][wb=0, wf1, wf2] - G[(i2, i1)][wf1] * G[(i3, i4)][wf2]

    if freqs is None: box frequency sampling
    else: sparse sampling

    """

    # g_ij[(i, j)] = data[iw]
    g_ij = {}
    for isp, sp in enumerate(spin_names):
        # data[iw, orb1, orb2]
        norb = gimp[sp].data.shape[1]
        assert norb == gimp[sp].data.shape[2]
        for o1, o2 in product(list(range(norb)), repeat=2):
            i = o1 + isp*norb
            j = o2 + isp*norb
            g_ij[(i, j)] = gimp[sp].data[:, o1, o2]

    assert g_ij[(0, 0)].shape[0] % 2 == 0
    w0 = g_ij[(0, 0)].shape[0] // 2

    def g(_i, _j, _w):
        return g_ij[(_i, _j)][w0 + _w]

    for (i1, i2, i3, i4), data in list(xloc.items()):
        # skip if g_ij has no data
        if (not (i2, i1) in g_ij) or (not (i3, i4) in g_ij):
            continue

        if freqs is None:
            # box sampling
            n_wf = data.shape[1]
            assert n_wf == data.shape[2]
            assert n_wf % 2 == 0

            for if1, if2 in product(list(range(n_wf)), repeat=2):
                wf1 = if1 - n_wf//2
                wf2 = if2 - n_wf//2
                data[0, if1, if2] -= g(i2, i1, wf1) * g(i3, i4, wf2)
        else:
            # sparse sampling
            assert freqs.shape[0] == data.shape[0]
            for j, (wb, wf1, wf2) in enumerate(freqs):
                if wb==0:
                    data[j] -= g(i2, i1, wf1) * g(i3, i4, wf2)


def gen_sparse_freqs_fix_boson(boson_freq, Lambda, sv_cutoff):
    """
    Generate sparse frequency points (wf1, wf2) for fixed bosonic freq wb

    Parameters
    ----------
    boson_freqs: (int) bosonic Matsubara frequency
    Lambda: (float)
    sv_cutoff: (float)

    Returns
    -------
    freqs: list of tuple [(wf1, wf2),]

    """
    # import irbasis
    from packaging.version import parse
    import irbasis_util
    if parse(irbasis_util.__version__) < parse('0.8'):
        raise RuntimeError('Please update irbasis-utility library!')
    from irbasis_util.four_point_ph_view import FourPointPHView

    print("  sparse_freqs: boson_freq = %3d," % boson_freq, end="")
    beta_dummy = 1.0
    basis = FourPointPHView(boson_freq, Lambda, beta_dummy, sv_cutoff)
    sp = basis.sampling_points_matsubara(whichl=basis.Nl-1)

    print("  # of points = %d" % len(sp))
    return sp

def gen_sparse_freqs_3D(Lambda, sv_cutoff):
    """
    Generate sparse frequency points (wb, wf1, wf2)

    Parameters
    ----------
    Lambda: (float)
    sv_cutoff: (float)

    Returns
    -------
    freqs: list of tuple [(wb, wf1, wf2),]

    """
    # import irbasis
    from irbasis_util.four_point import FourPoint, to_PH_convention

    beta_dummy = 1.0
    basis = FourPoint(Lambda, beta_dummy, sv_cutoff)
    # sp_ffff are in the convetion of four fermion frequencies
    sp_ffff = basis.sampling_points_matsubara(whichl=basis.Nl-1)

    print("  # of points = %d" % len(sp_ffff))
    return [to_PH_convention(s) for s in sp_ffff]


def gen_sparse_freqs(boson_freqs, Lambda, sv_cutoff):
    """
    Generate sparse frequency points (wb, wf1, wf2)

    Parameters
    ----------
    boson_freqs: (int) Number of bosonic Matsubara frequencies
    Lambda: (float)
    sv_cutoff: (float)

    Returns
    -------
    freqs: numpy.array of shape=(N, 3)

    """

    # [wb][i] = tuple(wf1, wf2)
    freqs2_list = [gen_sparse_freqs_fix_boson(boson_freq, Lambda, sv_cutoff)
                   for boson_freq in range(boson_freqs)]

    # convert to numpy.array of size [N, 3]
    freqs = []
    for wb, freqs2 in enumerate(freqs2_list):
        freqs3 = [(wb, wf1, wf2) for wf1, wf2 in freqs2]
        freqs.extend(freqs3)
    print(" Total sampling points = %d" % len(freqs))

    return numpy.array(freqs)


class SaveBSE:

    def __init__(self, h5_file, bse_info, n_corr_shells, n_flavors, use_spin_orbit, nonlocal_order_parameter, beta,
                 spin_names, bse_grp='', sparse_grp='bse_sparse'):

        # store all arguments for dump
        self.args = copy.deepcopy(locals())
        del self.args['self']

        from dft_tools.index_pair import IndexPair, IndexPair2

        self.spin_names = spin_names
        self.use_spin_orbit = use_spin_orbit
        self.n_flavors = n_flavors

        only_diagonal = not nonlocal_order_parameter

        n_block = 2 if not use_spin_orbit else 1
        n_inner = n_flavors // n_block
        self.n_orb = n_flavors // 2
        assert n_block == len(spin_names)

        # NOTE: change the order of spins in HDF5 to meet SumkDFTChi
        self.block2 = IndexPair2(list(range(n_corr_shells)), sorted(spin_names), only_diagonal1=only_diagonal)
        self.inner2 = IndexPair(list(range(n_inner)), convert_to_int=True)
        print(" block2 namelist =", self.block2.namelist)
        print(" inner2 namelist =", self.inner2.namelist)

        self.h5_file = h5_file
        self.bse_grp = bse_grp
        h5bse = h5BSE(self.h5_file, self.bse_grp)
        if bse_info == 'check':
            # check equivalence of info
            #assert compare_str_list(h5bse.get(key=('block_name',)), self.block2.namelist)
            #assert compare_str_list(h5bse.get(key=('inner_name',)), self.inner2.namelist)
            assert h5bse.get(key=('beta',)) == beta
        elif bse_info == 'save':
            # save info
            h5bse.save(key=('block_name',), data=self.block2.namelist)
            h5bse.save(key=('inner_name',), data=self.inner2.namelist)
            h5bse.save(key=('beta',), data=beta)
        else:
            raise ValueError("bse_info =", bse_info)

    def _save_common(self, xloc_ijkl, icrsh, key_type):

        h5bse = h5BSE(self.h5_file, self.bse_grp)

        n_inner2 = len(self.inner2.namelist)
        n_w2b = xloc_ijkl[list(xloc_ijkl.keys())[0]].shape[0]

        def decompose_index(index):
            spn = index // self.n_orb
            orb = index % self.n_orb
            return self.spin_names[spn], orb

        # read X_loc data and save into h5 file
        for wb in range(n_w2b):
            # boson freq
            xloc_bse = {}
            for (i1, i2, i3, i4), data in list(xloc_ijkl.items()):
                # (wb, wf1, wf2) --> (wf1, wf2)
                data_wb = data[wb]

                if not self.use_spin_orbit:
                    s1, o1 = decompose_index(i1)
                    s2, o2 = decompose_index(i2)
                    s3, o3 = decompose_index(i3)
                    s4, o4 = decompose_index(i4)
                else:
                    s1, o1 = 0, i1
                    s2, o2 = 0, i2
                    s3, o3 = 0, i3
                    s4, o4 = 0, i4

                s12 = self.block2.get_index(icrsh, s1, icrsh, s2)
                s34 = self.block2.get_index(icrsh, s3, icrsh, s4)
                inner12 = self.inner2.get_index(o1, o2)
                inner34 = self.inner2.get_index(o3, o4)

                if (s12, s34) not in xloc_bse:
                    xloc_bse[(s12, s34)] = numpy.zeros((n_inner2, n_inner2) + data_wb.shape, dtype=complex)
                xloc_bse[(s12, s34)][inner12, inner34] = data_wb  # copy

            # save
            h5bse.save(key=(key_type, wb), data=xloc_bse)

    def save_xloc(self, xloc_ijkl, icrsh):
        self._save_common(xloc_ijkl, icrsh, 'X_loc')

    def save_chiloc(self, chiloc_ijkl, icrsh):
        self._save_common(chiloc_ijkl, icrsh, 'chi_loc')

    def save_gamma0(self, u_mat, icrsh):

        assert u_mat.shape == (self.n_flavors, )*4

        h5bse = h5BSE(self.h5_file, self.bse_grp)

        # transpose U matrix into p-h channel, c_1^+ c_2 c_4^+ c_3
        u_mat_ph1 = u_mat.transpose(0, 2, 3, 1)
        u_mat_ph2 = u_mat.transpose(0, 3, 2, 1)

        if not self.use_spin_orbit:
            u_mat_ph1 = u_mat_ph1.reshape((2, self.n_orb)*4)
            u_mat_ph2 = u_mat_ph2.reshape((2, self.n_orb)*4)

            gamma0 = {}
            for s1, s2, s3, s4 in product(list(range(2)), repeat=4):
                gamma0_orb = - u_mat_ph1[s1, :, s2, :, s3, :, s4, :] + u_mat_ph2[s1, :, s2, :, s3, :, s4, :]

                # skip if zero
                if numpy.linalg.norm(gamma0_orb) == 0:
                    continue

                s12 = self.block2.get_index(icrsh, self.spin_names[s1], icrsh, self.spin_names[s2])
                s34 = self.block2.get_index(icrsh, self.spin_names[s3], icrsh, self.spin_names[s4])

                gamma0[(s12, s34)] = gamma0_orb.reshape((self.n_orb**2, )*2)

            h5bse.save(key=('gamma0', ), data=gamma0)
        else:
            gamma0_inner = - u_mat_ph1 + u_mat_ph2
            gamma0_inner = gamma0_inner.reshape((len(self.inner2.namelist), )*2)
            block_index = self.block2.get_index(icrsh, self.spin_names[0], icrsh, self.spin_names[0])
            h5bse.save(key=('gamma0', ), data={(block_index, block_index): gamma0_inner})

    def save_sparse_info(self, freqs):
        grp = self.args['sparse_grp']
        with HDFArchive(self.args['h5_file'], 'a') as f:
            if grp not in f:
                f.create_group(grp)
            # save all arguments in __init__ to re-instantiate
            f[grp]['args'] = self.args
            # save list of frequencies (sampling points)
            f[grp]['freqs'] = freqs

    def save_xloc_sparse(self, xloc_ijkl, icrsh):
        grp = self.args['sparse_grp']
        with HDFArchive(self.args['h5_file'], 'a') as f:
            f[grp][str(icrsh)] = xloc_ijkl

    @staticmethod
    def get_sparse_info(h5_file, sparse_grp='bse_sparse'):
        with HDFArchive(h5_file, 'r') as f:
            # if sparse_grp not in f:
            # save all arguments in __init__ to re-instantiate
            return f[sparse_grp]['args']



class DMFTBSESolver(DMFTCoreSolver):
    def __init__(self, seedname, params, output_file='', output_group='dmft_out'):
        super(DMFTBSESolver, self).__init__(seedname, params, output_file, output_group, read_only=True, restart=True)

    def _calc_bse_x0q(self):
        """
        Calc X_0(q)
        """
        print("\n--- dcore_bse - X_0(q)")
        if self._params['bse']['skip_X0q_if_exists'] and os.path.exists(self._params['bse']['h5_output_file']):
            print(" skip")
            return

        from .lattice_models import create_lattice_model

        lattice_model = create_lattice_model(self._params)

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'bse'
        params['mu'] = self._chemical_potential
        params['list_wb'] = numpy.arange(self._params['bse']['num_wb']).tolist()
        params['n_wf_G2'] = self._params['bse']['num_wf']
        params['div'] = lattice_model.nkdiv()
        params['bse_h5_out_file'] = os.path.abspath(self._params['bse']['h5_output_file'])
        params['use_temp_file'] = self._params['bse']['use_temp_file']
        if self._params['bse']['X0q_qpoints_saved'] == 'quadrant':
            params['X0q_qpoints_saved'] = 'quadrant'
        else:
            q_points = []
            with open(self._params['bse']['X0q_qpoints_saved'], 'r') as f:
                for line in f:
                    q_str = line.split()[1]
                    q_points.append(tuple(map(int, q_str.split('.'))))
            params['X0q_qpoints_saved'] = q_points
        r = run_sumkdft(
            'SumkDFTWorkerBSE',
            os.path.abspath(self._seedname+'.h5'), './work/sumkdft_bse', self._mpirun_command, params)

    def _calc_bse_xloc(self):
        """
        Calc X_loc and U matrix
        """

        # NOTE:
        #   n_flavors may depend on ish, but present BSE code does not support it
        def calc_num_flavors(_ish):
            return numpy.sum([len(indices) for indices in list(self._gf_struct[_ish].values())])
        n_flavors = calc_num_flavors(0)
        for ish in range(self._n_inequiv_shells):
            assert n_flavors == calc_num_flavors(ish)

        #
        # init for saving data into HDF5
        #
        print("\n--- dcore_bse - invoking h5BSE...")
        bse = SaveBSE(
            n_corr_shells=self._n_corr_shells,
            h5_file=os.path.abspath(self._params['bse']['h5_output_file']),
            bse_info='check',
            nonlocal_order_parameter=False,
            use_spin_orbit=self._use_spin_orbit,
            beta=self._beta,
            n_flavors=n_flavors,
            spin_names=self.spin_block_names)

        #
        # save U matrix for RPA
        #
        print("\n--- dcore_bse - U matrix")
        for ish in range(self._n_inequiv_shells):
            bse.save_gamma0(self._Umat[ish], icrsh=self._sk.inequiv_to_corr[ish])

        # FIXME:
        #     Saving data should be done for all **correlated_shells** (not for inequiv_shells)
        #     Namely, we need a loop for correlated shells when n_inequiv_shells < n_corr_shells
        #assert self._n_inequiv_shells == self._n_corr_shells

        #
        # X_loc
        #
        print("\n--- dcore_bse - X_loc")
        if self._params['bse']['skip_Xloc']:
            print(" skip")
            return

        Gloc_iw_sh, _ = self.calc_Gloc()
        solver_name = self._params['impurity_solver']['name']

        # generate sampling points of Matsubara frequencies
        print("\nFrequency sampling: box")
        freqs = None

        for ish in range(self._n_inequiv_shells):
            print("\nSolving impurity model for inequivalent shell " + str(ish) + " ...")
            sys.stdout.flush()
            x_loc, chi_loc, g_imp = calc_g2_in_impurity_model(solver_name, self._solver_params, self._mpirun_command,
                                                              self._params["impurity_solver"]["basis_rotation"],
                                                              self._Umat[ish], self._gf_struct[ish],
                                                              self._beta, self._n_iw,
                                                              self._sh_quant[ish].Sigma_iw, Gloc_iw_sh[ish],
                                                              self._params['bse']['num_wb'],
                                                              self._params['bse']['num_wf'], ish, freqs=freqs)

            subtract_disconnected(x_loc, g_imp, self.spin_block_names, freqs=freqs)

            # save X_loc, chi_loc
            for icrsh in range(self._n_corr_shells):
                if ish == self._sk.corr_to_inequiv[icrsh]:
                    # X_loc
                    bse.save_xloc(x_loc, icrsh=icrsh)
                    # chi_loc
                    if chi_loc is not None:
                        bse.save_chiloc(chi_loc, icrsh=icrsh)

    def calc_bse(self):
        """
        Compute data for BSE
        """
        self._calc_bse_x0q()
        self._calc_bse_xloc()


    def fit_Xloc(self, save_only):
        """
        Fit sparse data of Xloc
        """

        # Prepare input files
        print('')
        print('Fitting Xloc...')
        h5_path = os.path.abspath(self._params['bse']['h5_output_file'])
        work_dir = './work/sparse_fit-D{}'.format(self._params['bse']['sparse_D'])

        if not save_only:
            make_empty_dir(work_dir)

        cwd_org = os.getcwd()
        os.chdir(work_dir)

        if not save_only:
            t1 = time.time()
            for b in range(self._params['bse']['num_wb']):
                print('  wb={}...'.format(b))
                sys.stdout.flush()
                with open('./output-wb{}.txt'.format(b), 'w') as fout:
                    commands = [sys.executable, "-m", "dcore.sparse_sampling.ph"]
                    commands.extend(['--Lambda', self._params['bse']['sparse_Lambda']])
                    commands.extend(['--svcutoff', self._params['bse']['sparse_sv_cutoff']])
                    commands.extend(['--D', self._params['bse']['sparse_D']])
                    commands.extend(['--niter', self._params['bse']['sparse_niter']])
                    commands.extend(['--bfreq', b])
                    commands.extend(['--num_wf', self._params['bse']['num_wf']])
                    commands.extend(['--rtol', self._params['bse']['sparse_rtol']])
                    commands.extend(['--alpha', self._params['bse']['sparse_alpha']])
                    commands.append(h5_path)
                    commands = list(map(str, commands))
                    launch_mpi_subprocesses(self._mpirun_command, commands, fout)
                sys.stdout.flush()
            t2 = time.time()
            print('Fit ran for {} seconds'.format(t2-t1))
            sys.stdout.flush()

        # save interpolated data for BSE
        print('\n saving Xloc for BSE...')
        bse_args = SaveBSE.get_sparse_info(h5_path)
        bse = SaveBSE(**bse_args)

        for ish in range(self._n_inequiv_shells):
            # print('  ish={}...'.format(ish))

            x_loc = {}
            num_wb = self._params['bse']['num_wb']
            num_wf = self._params['bse']['num_wf']

            # get X_loc
            with h5py.File(h5_path, 'r') as f:
                for b in range(num_wb):
                # print('  wb={}...'.format(b))
                    prefix = '/bse_sparse/interpolated/{}/wb{}/D{}'.format(ish, b, self._params['bse']['sparse_D'])
                    assert prefix in f
                    for key in list(f[prefix].keys()):
                        # print(key)
                        key_tuple = ast.literal_eval(key)  # convert string to tuple
                        assert isinstance(key_tuple, tuple)

                        if key_tuple not in x_loc:
                            x_loc[key_tuple] = numpy.zeros((num_wb, 2*num_wf, 2*num_wf), dtype=complex)

                        x_wb_fix = float_to_complex_array(f[prefix + '/' + key][()])
                        assert x_wb_fix.shape == (2*num_wf, 2*num_wf)
                        x_loc[key_tuple][b, :, :] = x_wb_fix

            # save X_loc
            for icrsh in range(self._n_corr_shells):
                if ish == self._sk.corr_to_inequiv[icrsh]:
                    bse.save_xloc(x_loc, icrsh=icrsh)

        os.chdir(cwd_org)


def dcore_bse(filename, np=1):
    """
    Main routine for the BSE post-processing tool

    Parameters
    ----------
    filename : string
        Input-file name
    """

    print("\n############  Reading Input File  #################\n")
    print("  Input File Name : ", filename)

    #
    # Construct a parser with default values
    #
    pars = create_parser()

    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    seedname = p["model"]["seedname"]
    p["mpi"]["num_processes"] = np

    #
    # Load DMFT data
    #
    solver = DMFTBSESolver(seedname, p, output_file=seedname + '.out.h5')
    if solver.iteration_number == 0:
        raise RuntimeError("Number of iterations is zero!")
    print("Number of iterations :", solver.iteration_number)

    #
    # Compute data for BSE
    #
    solver.calc_bse()


    #
    # Finish
    #
    print("\n#################  Done  #####################\n")


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version, print_header

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_bse.py',
        description='Post-processing script for dcore (bse).',
        usage='$ dcore_bse input --np 4',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description(),
        add_help=True)
    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name."
                        )
    parser.add_argument('--np', help='Number of MPI processes', required=True)
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)
    dcore_bse(args.path_input_file, int(args.np))
