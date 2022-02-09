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


# here, any module depending on triqs.mpi or mpi4py should not be imported.
import sys
import os
import re
import time
import numpy
import copy
import ast
import h5py
import builtins
from typing import List

from ._dispatcher import dyson, make_zero_tail
from .program_options import *
from .sumkdft_workers.launcher import run_sumkdft

from .sumkdft_compat import SumkDFTCompat

from .tools import *
from .dc import hf_dc

from . import impurity_solvers
from .symmetrizer import pm_symmetrizer

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# def _total_density(bgf):
#     """Compute total density of a BlockGf instance assuming 1/iw"""
#     assert isinstance(bgf, BlockGf)
#     total_density = 0
#     for _, g_iw in bgf:
#         km = make_zero_tail(g_iw, 2)
#         km[1] = numpy.eye(g_iw.target_shape[0])
#         try:
#             total_density += g_iw.total_density(km).real
#         except RuntimeError:
#             from dcorelib.triqs_compat.gf import GfImFreq
#             g_iw_ = GfImFreq.from_triqs(g_iw)
#             total_density += g_iw_.total_density().real
#     return total_density

def _total_density(bgf, tail_fit=True):
    assert isinstance(bgf, BlockGf)
    if tail_fit:
        return bgf.total_density().real
    else:
        return calc_total_density(bgf).real

def _density(bgf, tail_fit=True):
    assert isinstance(bgf, BlockGf)
    if tail_fit:
        return bgf.density()
    else:
        return calc_density_matrix(bgf)


def __gettype(name):
    t = getattr(builtins, name)
    if isinstance(t, type):
        return t
    raise ValueError(name)

def create_solver_params(ini_dict):
    """
    Parse a dict and create parameters for an impurity solver.
    In dict, keyname should be parameter_name{python_type_name} (e.g. max_time{int}).

    :param ini_dict: a dict object containing keys and values read from *.ini
    :return: a dict object containing parameters
    """
    solver_params = {}
    for k, v in list(ini_dict.items()):
        r = re.compile('^(.*)\{(.*)\}$')
        try:
            m = r.search(k)
            if m is None:
                continue

            param_name, param_type_str = m.group(1), m.group(2)
            param_type = __gettype(param_type_str)
        except RuntimeError:
            raise RuntimeError("Unknown type or unrecognized format : " + k)
        solver_params[param_name] = param_type(v)

    return solver_params


class ShellQuantity(object):
    """

    Stores local quantities defined at a correlated shell
    All Gf-like objects are initialized to zero.
    All quantities are represented in LOCAL COORDINATE SYSTEM.

    """

    def __init__(self, gf_struct, beta, n_iw):
        # Self-energy
        self.Sigma_iw = make_block_gf(GfImFreq, gf_struct, beta, n_iw)
        for name, g in self.Sigma_iw:
            g.zero()


def _load_symm_generators(input_str, spin_orbit, dim_sh):
    """
    Load generators for symmetrization

    Parameters
    ----------
    input_str: input string
    spin_orbit: bool
      If spin_orbit is on or not
    dim_sh: int array
      dim of shells

    Returns
    -------
    List of generators

    """
    nsh = len(dim_sh)
    results = [[]] * nsh

    if input_str == '' or input_str == 'None':
        return results

    try:
        generators_file = ast.literal_eval(input_str)
    except:
        raise RuntimeError('Error in parsing {}'.format(input_str))

    assert isinstance(generators_file, list)

    for gen in generators_file:
        assert isinstance(gen, tuple) or isinstance(gen, list)
        assert len(gen) == 2
        assert isinstance(gen[0], int)
        assert isinstance(gen[1], str)

        ish = gen[0]
        print('Reading symmetry generator for ish={} from {}'.format(ish, gen[1]))
        if spin_orbit:
            gen_mat = numpy.zeros((1, dim_sh[ish], dim_sh[ish]), dtype=complex)
            read_potential(gen[1], gen_mat)
            results[ish].append({'ud' : gen_mat[0,:,:]})
        else:
            gen_mat = numpy.zeros((2, dim_sh[ish], dim_sh[ish]), dtype=complex)
            read_potential(gen[1], gen_mat)
            results[ish].append({'up' : gen_mat[0,:,:], 'down' : gen_mat[1,:,:]})
    return results


def solve_impurity_model(solver_name, solver_params, mpirun_command, basis_rot, Umat, gf_struct, beta, n_iw, Sigma_iw, Gloc_iw, mesh, ish, work_dir):
    """

    Solve an impurity model

    If mesh is not None, Sigma_w will be computed. Otherwise, None will be returned as Sigma_w.

    """

    assert isinstance(basis_rot, str)

    Solver = impurity_solvers.solver_classes[solver_name]

    raise_if_mpi_imported()

    sol = Solver(beta, gf_struct, Umat, n_iw)

    if not mesh is None and not Solver.is_gf_realomega_available():
        raise RuntimeError("Error: requested real-frequency quantities for an imaginary-time solver.")

    # Correct?
    # G0_iw^{-1} = Gloc_iw + Sigma_iw
    G0_iw = dyson(Sigma_iw=Sigma_iw, G_iw=Gloc_iw)
    diff = make_hermite_conjugate(G0_iw)
    if diff > 1e-8:
        print('Warning G0_iw is not hermite conjugate: {}'.format(diff))
    sol.set_G0_iw(G0_iw)

    s_params = copy.deepcopy(solver_params)
    s_params['random_seed_offset'] = 1000 * ish

    work_dir_org = os.getcwd()
    make_empty_dir(work_dir)
    os.chdir(work_dir)

    if not mesh is None:
        s_params['calc_Sigma_w'] = True
        s_params['omega_min'], s_params['omega_max'], s_params['n_omega'] = mesh

    # Solve the model
    rot = impurity_solvers.compute_basis_rot(basis_rot, sol)
    sol.solve(rot, mpirun_command, s_params)

    os.chdir(work_dir_org)

    # Read & save local quantities
    # Change from DCore v1:
    #      Local impurity Green's function is saved as "Gimp_iw" in DCore v2.
    #      This is intended to distinguish the local impurity Green's function from the one computed by SumkDFT.

    r = {
        'Sigma_iw': sol.get_Sigma_iw(),
        'Gimp_iw': sol.get_Gimp_iw(),
        'Sigma_w': sol.get_Sigma_w(),
        'quant_to_save': sol.quant_to_save,
    }
    return r


def calc_dc(dc_type, u_mat, dens_mat, spin_block_names, use_spin_orbit):

    # dim_tot is the dimension of spin x orbit for SO = 1 or that of orbital for SO=0
    # dim_tot = self._dim_sh[ish]
    # TODO: num_orb can be deleted, if inverse transformation of .tools.to_spin_full_U_matrix is implemented}
    num_orb = int(u_mat.shape[0] / 2)

    dc_imp_sh = {}  # This will be returned
    #
    # Hartree-Fock contribution
    #
    if dc_type == "HF_DFT" or dc_type == "HF_imp":
        if use_spin_orbit:
            dim_tot = dens_mat["ud"].shape[0]  # 2 * num_orb
            dc_imp_sh["ud"] = numpy.zeros((dim_tot, dim_tot), numpy.complex_)
            dc_imp_sh["ud"] += numpy.einsum("ijkl, jl->ik", u_mat, dens_mat["ud"], optimize=True)  # Hartree
            dc_imp_sh["ud"] -= numpy.einsum("ijkl, jk->il", u_mat, dens_mat["ud"], optimize=True)  # Fock
        else:
            for sp1 in spin_block_names:
                u_mat_reduce = u_mat[0:num_orb, 0:num_orb, 0:num_orb, 0:num_orb]  # TODO: make a function
                dens_mat_tot = sum(dens_mat.values())  # spin-sum of density matrix

                dc_imp_sh[sp1] = numpy.zeros((num_orb, num_orb), numpy.complex_)
                dc_imp_sh[sp1] += numpy.einsum("ijkl, jl->ik", u_mat_reduce, dens_mat_tot, optimize=True)  # Hartree
                dc_imp_sh[sp1] -= numpy.einsum("ijkl, jk->il", u_mat_reduce, dens_mat[sp1], optimize=True)  # Fock
    #
    # Fully-Localized Limit (FLL)
    #
    elif dc_type == "FLL":
        if use_spin_orbit:
            raise NotImplementedError
        else:
            u_sum = numpy.einsum("ijij", u_mat[0:num_orb, 0:num_orb, 0:num_orb, 0:num_orb], optimize=True)
            j_sum = numpy.einsum("ijji", u_mat[0:num_orb, 0:num_orb, 0:num_orb, 0:num_orb], optimize=True)

            u_ave = u_sum / num_orb**2  # U
            u_j_ave = (u_sum - j_sum) / (num_orb*(num_orb-1))  # U-J
            j_ave = u_ave - u_j_ave  # J

            n_sp = {sp: numpy.sum(numpy.diag(dens_mat[sp])).real for sp in spin_block_names}
            n_tot = sum(n_sp.values())

            print("\n      FLL formula")
            print("        U_ave = {0:.3f}".format(u_ave))
            print("        J_ave = {0:.3f}".format(j_ave))
            print("        n_sp  = {}".format(n_sp))
            print("        n_tot = {}".format(n_tot))

            dc_imp_sh = {}
            for sp in spin_block_names:
                dc = u_ave * (n_tot - 0.5) - j_ave * (n_sp[sp] - 0.5)
                dc_imp_sh[sp] = numpy.identity(num_orb) * dc
    else:
        raise ValueError("Here should not be reached")

    return dc_imp_sh


class DMFTCoreSolver(object):
    def __init__(self, seedname, params, output_file='', output_group='dmft_out', read_only=False, restart=False):
        """

        if read_only is True, no further SCF loops are performed. Chemical potential is fixed to that in the HDF5 file.

        :param seedname:
        :param params:
        :param output_file:
        :param output_group:
        :param read_only:
        """
        # Set up
        self._seedname = seedname
        self._params = copy.deepcopy(params)
        self._output_file = seedname+'.out.h5' if output_file == '' else output_file
        self._output_group = output_group
        self._use_spin_orbit = False

        self._beta = float(params['system']['beta'])
        self._n_iw = int(params['system']['n_iw'])  # Number of Matsubara frequencies

        # MPI commands
        if 'mpi' in params and 'num_processes' in params['mpi']:
            self._mpirun_command = params['mpi']['command'].replace('#', str(params['mpi']['num_processes']))

        self._read_only = read_only
        if read_only:
            assert restart == True

        #
        # Read dft input data
        #
        self._sk = SumkDFTCompat(hdf_file=seedname+'.h5')
        sk = self._sk

        self._use_spin_orbit = sk.SO != 0

        # Number of inequivalent shell
        self._n_inequiv_shells = sk.n_inequiv_shells

        # Number of correlated shell
        self._n_corr_shells = sk.n_corr_shells

        # Dimension of an inequivalent shell or a correlated shell
        # 'dim' is the dimension of either 'ud', 'up' or 'down' sectors for a correlated shell
        # For 'ud', dim is the size of spin times orbital.
        self._dim_corr_sh = numpy.array([sk.corr_shells[icrsh]['dim'] for icrsh in range(self._n_corr_shells)])
        self._dim_sh = numpy.array([self._dim_corr_sh[sk.inequiv_to_corr[ish]] for ish in range(self._n_inequiv_shells)])

        if self._use_spin_orbit:
            self._spin_block_names = ['ud']
        else:
            self._spin_block_names = ['up', 'down']

        # Structure of local Green's function at inequivalent shells
        self._gf_struct = []
        for ish in range(self._n_inequiv_shells):
            self._gf_struct.append(dict([(spn, numpy.arange(self._dim_sh[ish])) for spn in self._spin_block_names]))

        with HDFArchive(seedname+'.h5', 'r') as h:
            self._Umat = h["DCore"]["Umat"]

            # LocalPotential: convert to data structure similar to Sigma_iw_sh
            local_pot = h["DCore"]["LocalPotential"]
            def array2dict(array):
                # [sp, orb1, orb1] -> {sp_name: [orb1, orb2]}
                return {sp: array[i] for i, sp in enumerate(self._spin_block_names)}
            self._local_potential = [array2dict(local_pot_sh) for local_pot_sh in local_pot]

        # local quantities at ineuivalent shells
        self._sh_quant = [ShellQuantity(self._gf_struct[ish], self._beta, self._n_iw) for ish in range(self._n_inequiv_shells)]

        # DC correction at correlated shells
        self._dc_imp = []
        for icrsh in range(self._n_corr_shells):
            dc = {}
            for spn in self._spin_block_names:
                dc[spn] = numpy.zeros((self._dim_corr_sh[icrsh], self._dim_corr_sh[icrsh]), dtype=complex)
            self._dc_imp.append(dc)
        self._dc_energ = 0.0

        #
        # Read or set up seedname.out.h5
        #
        self._loaded_initial_self_energy = False
        if restart:
            if os.path.exists(self._output_file):
                self._read_output_file__restart()
                assert self._previous_runs >= 1
            else:
                print("Creating {}...".format(self._output_file))
                self._prepare_output_file__from_scratch()
                assert self._previous_runs == 0
        else:
            if os.path.exists(self._output_file):
                import shutil
                print("Moving {} to {}...".format(self._output_file, self._output_file + '.bak'))
                shutil.move(self._output_file, self._output_file + '.bak')

            self._prepare_output_file__from_scratch()
            assert self._previous_runs == 0

        if 'impurity_solver' in self._params:
            self._solver_params = create_solver_params(self._params['impurity_solver'])

        # physical quantities to save
        self._quant_to_save_latest = {}  # only the latest results are saved (overwritten)
        self._quant_to_save_history = {}  # histories are retained

        if not self._read_only:
            if self._params["control"]["time_reversal"]:
                norb_sh = self._dim_corr_sh//2 if self.use_spin_orbit else self._dim_corr_sh
                self._spin_symm = [
                    pm_symmetrizer(norb, self.use_spin_orbit, self._params["control"]["time_reversal_transverse"]) for norb in norb_sh]
            else:
                self._spin_symm = None

        self._sanity_check()

    def _read_output_file__restart(self):
        """
        Read data from & set up an output HDF5 file.
        """

        output_file, output_group = self._output_file, self._output_group

        # Set up a HDF file for output
        with HDFArchive(output_file, 'r') as f:
            ar = f[output_group]
            if 'iterations' not in ar:
                raise RuntimeError("Failed to restart the previous simulation! Data not found!")

            self._previous_runs = ar['iterations']
            if ar['iterations'] <= 0:
                raise RuntimeError("No previous runs to be loaded from " + output_file + "!")

            iterations = ar['iterations']

            print("Loading dc_imp and dc_energ... ")
            self._dc_imp = ar['dc_imp'][str(iterations)]
            self._dc_energ = ar['dc_energ'][str(iterations)]
            self._chemical_potential = ar['chemical_potential'][str(iterations)]
            if self._params['system']['fix_mu'] and self._chemical_potential != self._params['system']['mu']:
                raise RuntimeError('Wrong chemical potential {} found in {}! Given mu in system block = {}.'.format(self._chemical_potential, output_file, self._params['system']['mu']))

        # Load self-energy
        print("Loading Sigma_iw... ")
        with h5py.File(self._output_file, 'r') as ar:
            for ish in range(self._n_inequiv_shells):
                for bname, g in self._sh_quant[ish].Sigma_iw:
                    path = output_group + '/Sigma_iw/ite{}/sh{}/{}'.format(iterations, ish, bname)
                    load_giw(ar, path, g)
                diff = make_hermite_conjugate(self._sh_quant[ish].Sigma_iw)
                if diff > 1e-8:
                    print('Sigma_iw at shell {} is not exactly hermite: diff = {}, symmetrizing...'.format(ish, diff))

    def _prepare_output_file__from_scratch(self):
        """
        Set up an output HDF5 file.
        """

        assert not os.path.exists(self._output_file)

        self._chemical_potential = self._params['system']['mu']

        output_file, output_group = self._output_file, self._output_group

        # Set up a HDF file for output
        with HDFArchive(output_file, 'a') as f:
            if output_group in f:
                del f[output_group]

            f.create_group(output_group)
            f[output_group]['parameters'] = self._params

            #
            # Sub group for something
            #
            for gname in ['Sigma_iw', 'chemical_potential', 'dc_imp', 'dc_energ']:
                if not (gname in f[output_group]):
                    f[output_group].create_group(gname)

        self._previous_runs = 0

        # Set double-counting correction
        #     First, compute G_loc without self-energy
        if self._params['system']['with_dc']:
            print("@@@@@@@@@@@@@@@@@@@@@@@@  Double-Counting Correction  @@@@@@@@@@@@@@@@@@@@@@@@")
            _, dm0_sh = self.calc_G0loc()
            self.set_dc_imp(dm0_sh)

            # Add dc_imp to Sigma to cancel dc_imp in the first iteration
            for ish in range(self.n_inequiv_shells):
                for sp in self._spin_block_names:
                    self._sh_quant[ish].Sigma_iw[sp] << self._dc_imp[self._sk.inequiv_to_corr[ish]][sp]

        # Set initial value to self-energy
        if self._params["control"]["initial_static_self_energy"] != "None":
            print("@@@@@@@@@@@@@@@@@@@@@@@@  Setting initial value to self-energy @@@@@@@@@@@@@@@@@@@@@@@@")
            self._loaded_initial_self_energy = True
            init_se = set_potential(self._params["control"]["initial_static_self_energy"],
                                    "initial_static_self_energy",
                                    self.n_inequiv_shells, self._dim_sh, self.use_spin_orbit)

            for ish in range(self.n_inequiv_shells):
                for isp, sp in enumerate(self._spin_block_names):
                    # self._sh_quant[ish].Sigma_iw[sp] << init_se[ish][isp]
                    self._sh_quant[ish].Sigma_iw[sp].data[...] += init_se[ish][isp]

        elif self._params["control"]["initial_self_energy"] != "None":
            self._loaded_initial_self_energy = True
            Sigma_iw_sh_tmp = [self._sh_quant[ish].Sigma_iw for ish in range(self.n_inequiv_shells)]
            load_Sigma_iw_sh_txt(self._params["control"]["initial_self_energy"], Sigma_iw_sh_tmp, self.spin_block_names)
            print('')
            print('Loading {} ...'.format(self._params["control"]["initial_self_energy"]), end=' ')
            for ish in range(self.n_inequiv_shells):
                self._sh_quant[ish].Sigma_iw << Sigma_iw_sh_tmp[ish]
            print('Done')

        for ish in range(self.n_inequiv_shells):
            diff = make_hermite_conjugate(self._sh_quant[ish].Sigma_iw)
            if diff > 1e-8:
                print('Sigma(iwn) at ish = {} is not hermite (diff={}) but symmetrized.'.format(ish, diff))


    def _sanity_check(self):
        """
        Sanity checks
        """

        #if not numpy.allclose(self._sk.corr_to_inequiv, self._sk.inequiv_to_corr):
            #raise RuntimeError("corr_to_inequiv must be equal to inequiv_to_corr!")

        raise_if_mpi_imported()

        for icrsh in range(self._n_corr_shells):
            for sp1 in self._spin_block_names:
                if not numpy.allclose(self._dc_imp[icrsh][sp1], self._dc_imp[icrsh][sp1].conjugate().transpose()):
                    raise RuntimeError("dc_imp is not hermite!")

        for ish in range(self._n_inequiv_shells):
            if make_hermite_conjugate(self._sh_quant[ish].Sigma_iw, check_only=True) > 1e-8:
                raise RuntimeError("Sigma_iw is not hermite conjugate!")

    def _make_sumkdft_params(self):
        return {
            'beta'          : self._params['system']['beta'],
            'prec_mu'       : self._params['system']['prec_mu'],
            'with_dc'       : self._params['system']['with_dc'],
            'Sigma_iw_sh'   : [s.Sigma_iw for s in self._sh_quant],
            'potential'     : self._local_potential,
            'dc_imp'        : self._dc_imp,
            'dc_energ'      : self._dc_energ,
            'mu'            : self._chemical_potential,
            'adjust_mu'     : False,
            'no_tail_fit'   : self._params['system']['no_tail_fit'],
        }

    def calc_G0loc(self):
        """
        Compute the non-interacting lattice/local Green's function using SumkDFT.

        Although the chemical potential is adjusted on the fly so that the system has the correct number of electrons,
        self._chemical_potential is not updated.

        Return a list of Gloc_iw and density matrices for inequivalent shells
        """

        params = self._make_sumkdft_params()
        params['adjust_mu'] = True
        params['with_dc'] = False

        r = run_sumkdft(
            'SumkDFTWorkerGloc',
            os.path.abspath(self._seedname+'.h5'), './work/sumkdft_G0loc', self._mpirun_command, params)

        # Make sure Gloc_iw is hermite conjugate (for triqs 2.x)
        for ish, g in enumerate(r['Gloc_iw_sh']):
            diff = make_hermite_conjugate(g)
            if diff > 1e-8:
                print('Warning Gloc_iw at ish {} is not hermite conjugate: {}'.format(ish, diff))
        return r['Gloc_iw_sh'], r['dm_sh']


    def calc_Gloc(self):
        """
        Compute the lattice/local Green's function using SumkDFT.

        Return a list of Gloc_iw and density matrices for inequivalent shells
        """

        mu_old = self._chemical_potential

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'Gloc'
        if (not self._params['system']['fix_mu']) and (not self._read_only):
            params['adjust_mu'] = True
        r = run_sumkdft(
            'SumkDFTWorkerGloc',
            os.path.abspath(self._seedname+'.h5'), './work/sumkdft_Gloc', self._mpirun_command, params)

        if params['adjust_mu']:
            self._chemical_potential = r['mu']

        # Sanity check
        if self._params['system']['fix_mu'] or self._read_only:
            assert self._chemical_potential == mu_old

        # Make sure Gloc_iw is hermite conjugate (for triqs 2.x)
        for ish, g in enumerate(r['Gloc_iw_sh']):
            diff = make_hermite_conjugate(g)
            if diff > 1e-8:
                print('Warning Gloc_iw at ish {} is not hermite conjugate: {}'.format(ish, diff))
        return r['Gloc_iw_sh'], r['dm_sh']


    def print_density_matrix(self, dm_sh, smoment_sh):
        # dm_sh and smoment_sh are defined for inequivalent shells
        assert len(dm_sh) == len(smoment_sh) == self._n_inequiv_shells

        print("\nDensity Matrix")
        for ish in range(self._n_inequiv_shells):
            print("\n  Inequivalent Shell ", ish)
            for sp in self._spin_block_names:
                print("\n    Spin ", sp)
                for i1 in range(self._dim_sh[ish]):
                    print("          ", end="")
                    for i2 in range(self._dim_sh[ish]):
                        print("{0:.3f} ".format(dm_sh[ish][sp][i1, i2]), end="")
                    print("")
                evals, evecs = numpy.linalg.eigh(dm_sh[ish][sp])
                print('    Eigenvalues: ', evals)
            print('')
            print('    Magnetic moment (only spin contribution, S=1/2 gives 0.5)')
            print('      mx,my,mz= {} {} {}'.format(smoment_sh[ish][0], smoment_sh[ish][1], smoment_sh[ish][2]))


    def solve_impurity_models(self, Gloc_iw_sh, iteration_number, mesh=None):
        """

        Solve impurity models for all inequivalent shells

        :param Gloc_iw_sh:
        :param mesh: (float, float, int)
            (om_min, om_max, n_om)
        :return:
        """

        self._sanity_check()

        solver_name = self._params['impurity_solver']['name']
        Sigma_iw_sh = []
        Gimp_iw_sh = []
        Sigma_w_sh = []
        quant_to_save_sh = []
        for ish in range(self._n_inequiv_shells):
            print('')
            work_dir = 'work/imp_shell'+str(ish)+'_ite'+str(iteration_number)
            print('Solving impurity model for inequivalent shell {} in {}...'.format(ish, work_dir))
            print('')
            sys.stdout.flush()
            r = solve_impurity_model(solver_name, self._solver_params, self._mpirun_command,
                                     self._params["impurity_solver"]["basis_rotation"],
                                     self._Umat[ish], self._gf_struct[ish], self._beta, self._n_iw,
                                     self._sh_quant[ish].Sigma_iw, Gloc_iw_sh[ish], mesh, ish, work_dir)
            Sigma_iw = r['Sigma_iw']
            Gimp_iw = r['Gimp_iw']
            Sigma_w = r['Sigma_w']
            if make_hermite_conjugate(Sigma_iw) > 1e-8:
                print("Warning: Sigma_iw is not hermite conjugate!")
            if make_hermite_conjugate(Gimp_iw) > 1e-8:
                print("Warning: Gimp_iw is not hermite conjugate!")

            Sigma_iw_sh.append(Sigma_iw)
            if not mesh is None:
                Sigma_w_sh.append(Sigma_w)
            Gimp_iw_sh.append(Gimp_iw)

            quant_to_save_sh.append(r['quant_to_save'])

        # convert data structure of quant_to_save_sh from [sh][key] to [key][sh]
        # and add it to self._quant_to_save
        # TODO: maybe, for better visibility, return quant_to_save_sh and do the process below in do_steps()
        quant_to_save_dict = {key: [] for key in list(quant_to_save_sh[0].keys())}
        for quant_dict in quant_to_save_sh:
            for key, quant in list(quant_dict.items()):
                quant_to_save_dict[key].append(quant)
        self._quant_to_save_latest.update(quant_to_save_dict)

        # FIXME: DON'T CHANGE THE NUMBER OF RETURNED VALUES. USE DICT INSTEAD.
        if mesh is None:
            return Sigma_iw_sh, Gimp_iw_sh
        else:
            return Sigma_iw_sh, Gimp_iw_sh, Sigma_w_sh

    def set_dc_imp(self, dm_sh):
        """

        Compute Double-counting term (Hartree-Fock term)
        FIXME: PULL THIS OUT OF THIS CLASS

        """

        dc_type = self._params['system']['dc_type']
        print("\n  set DC (dc_type = {})".format(dc_type))

        def print_matrix(mat):
            dim = mat.shape[0]
            for i1 in range(dim):
                print("          ", end="")
                for i2 in range(dim):
                    print("{0:.3f} ".format(mat[i1, i2]), end="")
                print("")

        # dm_sh is defined for inequivalent shells
        assert len(dm_sh) == self._n_inequiv_shells

        # Loop over inequivalent shells
        _dc_imp = []
        for ish in range(self._n_inequiv_shells):
            u_mat = self._Umat[ish]
            dens_mat = dm_sh[ish]

            num_orb = int(u_mat.shape[0] / 2)
            # TODO: num_orb can be deleted, if inverse transformation of .tools.to_spin_full_U_matrix is implemented}

            # print U matrix, J matrix, density matrix
            print("\n    DC for inequivalent shell {0}".format(ish))
            print("\n      2-index U:".format(ish))
            print_matrix(numpy.einsum("ijij->ij", u_mat[0:num_orb, 0:num_orb, 0:num_orb, 0:num_orb]))
            print("\n      2-index J:".format(ish))
            print_matrix(numpy.einsum("ijji->ij", u_mat[0:num_orb, 0:num_orb, 0:num_orb, 0:num_orb]))
            print("\n      Local Density Matrix:".format(ish))
            for sp1 in self._spin_block_names:
                print("        Spin {0}".format(sp1))
                print_matrix(dens_mat[sp1][0:num_orb, 0:num_orb])

            # calculated DC
            _dc_imp.append(calc_dc(dc_type, u_mat, dens_mat, self.spin_block_names, self._use_spin_orbit))

            # print DC self energy
            print("\n      DC Self Energy:")
            for sp1 in self._spin_block_names:
                print("        Spin {0}".format(sp1))
                print_matrix(_dc_imp[ish][sp1])
            print("")

        # Now copy dc_imp to all correlated shells
        self._dc_imp = [_dc_imp[self._sk.corr_to_inequiv[icrsh]] for icrsh in range(self._n_corr_shells)]

    def do_steps(self, max_step):
        """

        Do more steps

        :param max_step: int
            Number of steps

        """

        assert not self._read_only

        previous_present = self._previous_runs > 0
        with_dc = self._params['system']['with_dc']
        dc_type = self._params['system']['dc_type']
        sigma_mix = self._params['control']['sigma_mix']  # Mixing factor of Sigma after solution of the AIM
        output_group = self._output_group

        symm_generators = _load_symm_generators(
            self._params['control']['symmetry_generators'],
            self._use_spin_orbit, self._dim_sh)

        # No tail fit in computing density if no_tail_fit=True
        tail_fit = not self._params['system']['no_tail_fit']

        def quantities_to_check():
            x = []
            # chemical potential
            x.append(self._chemical_potential)

            # renormalization factor
            w0 = self._n_iw  # number of positive Matsubara freqs = position of zeroth Matsubara freq
            # print(w0)
            for ish in range(self._n_inequiv_shells):
                for spn in self.spin_block_names:
                    for iorb in range(self._dim_sh[ish]):
                        # print(self._sh_quant[ish].Sigma_iw[spn].data.shape)
                        sigma0 = self._sh_quant[ish].Sigma_iw[spn].data[w0, iorb, iorb].imag
                        z = 1. / (1 - sigma0 / numpy.pi * self._beta)
                        x.append(z)
            return numpy.array(x, dtype=complex)

        x0 = quantities_to_check()

        converge_count = 0

        t0 = [time.time()] * 2
        def print_time(comment):
            t_now = time.time()
            print("\nWall Time : %.1f sec (%.1f sec for %s)" % (t_now - t0[0], t_now - t0[1], comment))
            t0[1] = t_now

        for iteration_number in range(self._previous_runs+1, self._previous_runs+max_step+1):
            self._sanity_check()

            print_time("start iteration %d" % iteration_number)
            sys.stdout.flush()
            print("")
            print("#####################################################################")
            print("########################  Iteration = %5d  ########################"%iteration_number)
            print("#####################################################################")
            print("")
            print("")
            print("@@@@@@@@@@@@@@@@@@@@@@@@  Chemical potential and G0_imp  @@@@@@@@@@@@@@@@@@@@@@@@")
            print("")
            sys.stdout.flush()

            # Compute Gloc_iw where the chemical potential is adjusted if needed
            Gloc_iw_sh, dm_sh = self.calc_Gloc()
            smoment_sh = spin_moments_sh(dm_sh)
            self.print_density_matrix(dm_sh, smoment_sh)
            self._quant_to_save_history['density_matrix'] = dm_sh
            self._quant_to_save_history['spin_moment'] = smoment_sh

            # Compute Total charge from G_loc
            charge_loc = [_total_density(Gloc_iw_sh[ish], tail_fit) for ish in range(self._n_inequiv_shells)]
            for ish, charge in enumerate(charge_loc):
                print("\n  Total charge of Gloc_{shell %d} : %.6f" % (ish, charge))
            self._quant_to_save_history['total_charge_loc'] = charge_loc

            print_time("calc_Gloc with chemical potential tuning")
            sys.stdout.flush()

            new_Sigma_iw, new_Gimp_iw = self.solve_impurity_models(Gloc_iw_sh, iteration_number)

            print_time("impurity problem")
            sys.stdout.flush()

            #
            # Solved. Now do post-processing:
            #

            # Compute Total charge from G_imp
            charge_imp = [_total_density(new_Gimp_iw[ish], tail_fit) for ish in range(self._n_inequiv_shells)]
            for ish, charge in enumerate(charge_imp):
                print("\n  Total charge of Gimp_{shell %d} : %.6f" % (ish, charge))
            self._quant_to_save_history['total_charge_imp'] = charge_imp

            # Symmetrize over spin components
            if self._spin_symm is not None:
                print("Averaging self-energy and impurity Green's function over spin components...")
                for ish in range(self._n_inequiv_shells):
                    new_Gimp_iw[ish] << self._spin_symm[ish](new_Gimp_iw[ish])
                    new_Sigma_iw[ish] << self._spin_symm[ish](new_Sigma_iw[ish])

            # Update Sigma_iw and Gimp_iw.
            # Mix Sigma if requested.
            if iteration_number > 1 or previous_present or self._loaded_initial_self_energy:
                for ish in range(self._n_inequiv_shells):
                    self._sh_quant[ish].Sigma_iw << sigma_mix * new_Sigma_iw[ish] \
                                + (1.0-sigma_mix) * self._sh_quant[ish].Sigma_iw
            else:
                for ish in range(self._n_inequiv_shells):
                    self._sh_quant[ish].Sigma_iw << new_Sigma_iw[ish]

            # Symmetrization
            for ish in range(self._n_inequiv_shells):
                if len(symm_generators[ish]) > 0:
                    self._sh_quant[ish].Sigma_iw << symmetrize(self._sh_quant[ish].Sigma_iw, symm_generators[ish])

            # update DC correction
            if with_dc and dc_type == "HF_imp":
                dm_imp = [_density(new_Gimp_iw[ish], tail_fit) for ish in range(self._n_inequiv_shells)]
                self.set_dc_imp(dm_imp)

            self._quant_to_save_history.update(chemical_potential=self._chemical_potential,
                                               dc_imp=self._dc_imp,
                                               dc_energ=self._dc_energ)

            # Save data to the hdf5 archive
            with HDFArchive(self._output_file, 'a') as ar:
                ar[output_group]['iterations'] = iteration_number

                # save histories
                for key, val in list(self._quant_to_save_history.items()):
                    if not key in ar[output_group]:
                        ar[output_group].create_group(key)
                    ar[output_group][key][str(iteration_number)] = val

                # save the latest results
                for key, val in list(self._quant_to_save_latest.items()):
                    ar[output_group][key] = val

            # Save the history of Sigma in DCore format
            with h5py.File(self._output_file, 'a') as ar:
                for ish in range(self._n_inequiv_shells):
                    for bname, g in self._sh_quant[ish].Sigma_iw:
                        path = output_group + '/Sigma_iw/ite{}/sh{}/{}'.format(iteration_number, ish, bname)
                        save_giw(ar, path, g)

            # Save Sigma in *.npz file
            self._save_sigma_iw(dm_sh)

            # convergence check
            tol = self._params["control"]["converge_tol"]
            if tol > 0:
                print("\nConvergence check")
                x1 = quantities_to_check()
                max_error = numpy.max(numpy.abs(x1 - x0))
                print(" | converge_tol = %.1e" %tol)
                print(" | max_error    = %.1e" %max_error)
                if max_error < tol:
                    converge_count += 1
                    print(" | convergence criterion satisfied. count={}".format(converge_count))
                else:
                    converge_count = 0
                    print(" | convergence criterion not satisfied. count={}".format(converge_count))
                x0 = x1

                if converge_count == self._params["control"]["n_converge"]:
                    print(" | converged --- iteration=%d" % iteration_number)
                    break

            sys.stdout.flush()

        self._previous_runs += max_step

    def _save_sigma_iw(self, dm_sh: List[numpy.ndarray]) -> None:
        """ Save Sigma(iw) for post processing/restart """
        data = {}
        idata = 0
        for ish in range(self._n_inequiv_shells):
            hf = hf_dc(dm_sh[ish], self._Umat[ish], self._use_spin_orbit)
            for bname, g in self._sh_quant[ish].Sigma_iw:
                data[f"data{idata}"] = g.data
                data[f"hartree_fock{idata}"] = hf[bname]
                idata += 1
        numpy.savez(self._seedname + "_sigma_iw.npz", **data)

    def chemical_potential(self, iteration_number):
        with HDFArchive(self._output_file, 'r') as ar:
            return ar[self._output_group]['chemical_potential'][str(iteration_number)]

    def density_matrix(self, iteration_number):
        with HDFArchive(self._output_file, 'r') as ar:
            return ar[self._output_group]['density_matrix'][str(iteration_number)]

    def spin_moment(self, iteration_number):
        with HDFArchive(self._output_file, 'r') as ar:
            return ar[self._output_group]['spin_moment'][str(iteration_number)]

    @property
    def n_inequiv_shells(self):
        return self._n_inequiv_shells

    @property
    def inequiv_to_corr(self):
        return self._sk.inequiv_to_corr

    @property
    def corr_to_inequiv(self):
        return self._sk.corr_to_inequiv

    @property
    def iteration_number(self):
        return self._previous_runs

    @property
    def spin_block_names(self):
        return self._spin_block_names

    @property
    def use_spin_orbit(self):
        return self._use_spin_orbit

    def inequiv_shell_info(self, ish):
        info = {}
        if self._use_spin_orbit:
            info['num_orb'] = self._dim_sh[ish]//2
        else:
            info['num_orb'] = self._dim_sh[ish]

        info['block_dim'] = self._dim_sh[ish]

        return info

    def corr_shell_info(self, ish):
        return self.inequiv_shell_info(self._sk.corr_to_inequiv[ish])

    def Sigma_iw_sh(self, iteration_number):
        Sigma_iw_sh = []
        with h5py.File(self._output_file, 'r') as ar:
            for ish in range(self._n_inequiv_shells):
                Sigma_iw = make_block_gf(GfImFreq, self._gf_struct[ish], self._beta, self._n_iw)
                for bname, g in Sigma_iw:
                    path = self._output_group + '/Sigma_iw/ite{}/sh{}/{}'.format(iteration_number, ish, bname)
                    load_giw(ar, path, g)
                Sigma_iw_sh.append(Sigma_iw)
        return Sigma_iw_sh


