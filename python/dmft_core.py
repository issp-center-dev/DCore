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
from __future__ import print_function

# here, any module depending on pytriqs.mpi or mpi4py should not be imported.
import sys
import os
import re
import time
import __builtin__
import numpy
import scipy
import copy
import ast
import h5py

from .program_options import *

from . import sumkdft

from tools import *

import impurity_solvers

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

def __gettype(name):
    t = getattr(__builtin__, name)
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
    for k, v in ini_dict.items():
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

    # Compute rotation matrix to the diagonal basis if supported
    if basis_rot == 'None':
        rot = None
    elif basis_rot == 'Hloc':
        rot = compute_diag_basis(G0_iw)
    else:
        if not os.path.exists(basis_rot):
            raise RuntimeError("Invalid basis_rot : {}".format(basis_rot))
        if sol.use_spin_orbit:
            rot = numpy.zeros((1, sol.n_flavors, sol.n_flavors), dtype=complex)
            read_potential(basis_rot, rot)
            rot = {'ud' : rot[0,:,:]}
        else:
            rot = numpy.zeros((2, sol.n_orb, sol.n_orb), dtype=complex)
            read_potential(basis_rot, rot)
            rot = {'up' : rot[0,:,:], 'down': rot[1,:,:]}
    s_params = copy.deepcopy(solver_params)
    s_params['random_seed_offset'] = 1000 * ish

    work_dir_org = os.getcwd()
    make_empty_dir(work_dir)
    os.chdir(work_dir)

    if not mesh is None:
        s_params['calc_Sigma_w'] = True
        s_params['omega_min'], s_params['omega_max'], s_params['n_omega'] = mesh

    # Solve the model
    sol.solve(rot, mpirun_command, s_params)

    os.chdir(work_dir_org)

    # Read & save local quantities
    # Change from DCore v1:
    #      Local impurity Green's function is saved as "Gimp_iw" in DCore v2.
    #      This is intended to distinguish the local impurity Green's function from the one computed by SumkDFT.

    return sol.get_Sigma_iw(), sol.get_Gimp_iw(), sol.get_Sigma_w()



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
        self._output_file = seedname+'.out.h5' if output_file is '' else output_file
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
        self._sk = sumkdft.SumkDFTCompat(hdf_file=seedname+'.h5')
        sk = self._sk

        self._use_spin_orbit = sk.SO != 0

        # Number of inequivalent shell
        self._n_inequiv_shells = sk.n_inequiv_shells

        # Number of correlated shell
        self._n_corr_shells = sk.n_corr_shells

        # Dimension of an inequivalent shell or a correlated shell
        # 'dim' is the dimension of either 'ud', 'up' or 'down' sectors for a correlated shell
        # For 'ud', dim is the size of spin times orbital.
        self._dim_corr_sh = [sk.corr_shells[icrsh]['dim'] for icrsh in range(self._n_corr_shells)]
        self._dim_sh = [self._dim_corr_sh[sk.inequiv_to_corr[ish]] for ish in range(self._n_inequiv_shells)]

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
            Gloc_iw_sh, dm_sh = self.calc_Gloc()
            self.set_dc_imp(dm_sh)

        # Set initial value to self-energy
        if self._params["control"]["initial_static_self_energy"] != "None":
            print("@@@@@@@@@@@@@@@@@@@@@@@@  Setting initial value to self-energy @@@@@@@@@@@@@@@@@@@@@@@@")
            init_se = set_potential(self._params["control"]["initial_static_self_energy"],
                                    "initial_static_self_energy",
                                    self.n_inequiv_shells, self._dim_sh, self.use_spin_orbit)

            for ish in range(self.n_inequiv_shells):
                for isp, sp in enumerate(self._spin_block_names):
                    self._sh_quant[ish].Sigma_iw[sp] << init_se[ish][isp]

        elif self._params["control"]["initial_self_energy"] != "None":
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

        for ish in range(self._n_inequiv_shells):
            for sp1 in self._spin_block_names:
                if not numpy.allclose(self._dc_imp[ish][sp1], self._dc_imp[ish][sp1].conjugate().transpose()):
                    raise RuntimeError("dc_imp is not hermite!")

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
        }

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
        r = sumkdft.run(os.path.abspath(self._seedname+'.h5'), './work/sumkdft', self._mpirun_command, params)

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


    def print_density_matrix(self, dm_sh):
        smoments = spin_moments_sh(dm_sh)
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
            print('      mx,my,mz= {} {} {}'.format(smoments[ish][0], smoments[ish][1], smoments[ish][2]))


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
        for ish in range(self._n_inequiv_shells):
            print('')
            work_dir = 'work/imp_shell'+str(ish)+'_ite'+str(iteration_number)
            print('Solving impurity model for inequivalent shell {} in {}...'.format(ish, work_dir))
            print('')
            sys.stdout.flush()
            Sigma_iw, Gimp_iw, Sigma_w = solve_impurity_model(solver_name, self._solver_params, self._mpirun_command,
                             self._params["impurity_solver"]["basis_rotation"], self._Umat[ish], self._gf_struct[ish],
                                 self._beta, self._n_iw,
                                 self._sh_quant[ish].Sigma_iw, Gloc_iw_sh[ish], mesh, ish, work_dir)
            if make_hermite_conjugate(Sigma_iw) > 1e-8:
                print("Warning: Sigma_iw is not hermite conjugate!")
            if make_hermite_conjugate(Gimp_iw) > 1e-8:
                print("Warning: Gimp_iw is not hermite conjugate!")

            Sigma_iw_sh.append(Sigma_iw)
            if not mesh is None:
                Sigma_w_sh.append(Sigma_w)
            Gimp_iw_sh.append(Gimp_iw)

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

        # Loop over inequivalent shells
        self._dc_imp = []
        for ish in range(self._n_inequiv_shells):
            u_mat = self._Umat[self._sk.inequiv_to_corr[ish]]

            # dim_tot is the dimension of spin x orbit for SO = 1 or that of orbital for SO=0
            dim_tot = self._dim_sh[ish]
            num_orb = int(u_mat.shape[0] / 2)

            dens_mat = dm_sh[self._sk.inequiv_to_corr[ish]]

            print("")
            print("    DC for inequivalent shell {0}".format(ish))
            print("\n      2-index U:".format(ish))
            for i1 in range(num_orb):
                print("          ", end="")
                for i2 in range(num_orb):
                    print("{0:.3f} ".format(u_mat[i1, i2, i1, i2]), end="")
                print("")
            print("\n      2-index J:".format(ish))
            for i1 in range(num_orb):
                print("          ", end="")
                for i2 in range(num_orb):
                    print("{0:.3f} ".format(u_mat[i1, i2, i2, i1]), end="")
                print("")

            print("\n      Local Density Matrix:".format(ish))
            for sp1 in self._spin_block_names:
                print("        Spin {0}".format(sp1))
                for i1 in range(num_orb):
                    print("          ", end="")
                    for i2 in range(num_orb):
                        print("{0:.3f} ".format(dens_mat[sp1][i1, i2]), end="")
                    print("")

            if self._use_spin_orbit:
                dc_imp_sh = {}
                dc_imp_sh["ud"] = numpy.zeros((dim_tot, dim_tot), numpy.complex_)
                for s1, i1, s2, i2 in product(range(2), range(num_orb), range(2), range(num_orb)):
                    #
                    # Hartree
                    #
                    dc_imp_sh["ud"][i1 + s1 * num_orb, i2 + s1 * num_orb] += numpy.sum(
                        u_mat[i1, 0:num_orb, i2, 0:num_orb] * dens_mat["ud"][s2 * num_orb:s2 * num_orb + num_orb,
                                                              s2 * num_orb:s2 * num_orb + num_orb]
                    )
                    #
                    # Exchange
                    #
                    dc_imp_sh["ud"][i1 + s1 * num_orb, i2 + s2 * num_orb] += numpy.sum(
                        u_mat[i1, 0:num_orb, 0:num_orb, i2]
                        * dens_mat["ud"][s2 * num_orb:s2 * num_orb + num_orb, s1 * num_orb:s1 * num_orb + num_orb]
                    )
                self._dc_imp.append(dc_imp_sh)
            else:
                dc_imp_sh = {}
                for sp1 in self._spin_block_names:
                    dc_imp_sh[sp1] = numpy.zeros((num_orb, num_orb), numpy.complex_)
                    for i1, i2 in product(range(num_orb), repeat=2):
                        #
                        # Hartree
                        #
                        for sp2 in self._spin_block_names:
                            dc_imp_sh[sp1][i1, i2] += \
                                numpy.sum(u_mat[i1, 0:num_orb, i2, 0:num_orb] * dens_mat[sp2][:, :])
                        #
                        # Exchange
                        #
                        dc_imp_sh[sp1][i1, i2] += \
                            - numpy.sum(u_mat[i1, 0:num_orb, 0:num_orb, i2] * dens_mat[sp1][:, :])
                self._dc_imp.append(dc_imp_sh)

            print("\n      DC Self Energy:")
            for sp1 in self._spin_block_names:
                print("        Spin {0}".format(sp1))
                for i1 in range(dim_tot):
                    print("          ", end="")
                    for i2 in range(dim_tot):
                        print("{0:.3f} ".format(self._dc_imp[ish][sp1][i1, i2]), end="")
                    print("")
            print("")


    def do_steps(self, max_step):
        """

        Do more steps

        :param max_step: int
            Number of steps

        """

        assert not self._read_only

        previous_present = self._previous_runs > 0
        with_dc = self._params['system']['with_dc']
        sigma_mix = self._params['control']['sigma_mix']  # Mixing factor of Sigma after solution of the AIM
        output_group = self._output_group

        symm_generators = _load_symm_generators(
            self._params['control']['symmetry_generators'],
            self._use_spin_orbit, self._dim_sh)

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

        t0 = time.time()
        for iteration_number in range(self._previous_runs+1, self._previous_runs+max_step+1):
            self._sanity_check()

            print("\nWall Time : %.1f sec" % (time.time() - t0))
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
            self.print_density_matrix(dm_sh)

            for ish in range(self._n_inequiv_shells):
                print("\n  Total charge of Gloc_{shell %d} : %.6f" % (ish, float(Gloc_iw_sh[ish].total_density())))

            print("\nWall Time : %.1f sec" % (time.time() - t0))
            sys.stdout.flush()

            new_Sigma_iw, new_Gimp_iw = self.solve_impurity_models(Gloc_iw_sh, iteration_number)

            print("\nWall Time : %.1f sec" % (time.time() - t0))
            sys.stdout.flush()

            # Solved. Now do post-processing:
            for ish in range(self._n_inequiv_shells):
                print("\nTotal charge of impurity problem : %.6f" % new_Gimp_iw[ish].total_density())

            # Symmetrize over spin components
            if self._params["control"]["time_reversal"]:
                print("Averaging self-energy and impurity Green's function over spin components...")

                if self._params["model"]["spin_orbit"]:
                    # TODO
                    raise Exception("Spin-symmetrization in the case with the spin-orbit coupling is not implemented")

                for ish in range(self._n_inequiv_shells):
                    symmetrize_spin(new_Gimp_iw[ish])
                    symmetrize_spin(new_Sigma_iw[ish])


            # Update Sigma_iw and Gimp_iw.
            # Mix Sigma if requested.
            if iteration_number > 1 or previous_present:
                for ish in range(self._n_inequiv_shells):
                    self._sh_quant[ish].Sigma_iw << sigma_mix * new_Sigma_iw[ish] \
                                + (1.0-sigma_mix) * self._sh_quant[ish].Sigma_iw
            else:
                for ish in range(self._n_inequiv_shells):
                    self._sh_quant[ish].Sigma_iw << new_Sigma_iw[ish]

            # Symmetrization
            for isn in range(self._n_inequiv_shells):
                if len(symm_generators[ish]) > 0:
                    self._sh_quant[ish].Sigma_iw << symmetrize(self._sh_quant[ish].Sigma_iw, symm_generators[ish])

            # Write data to the hdf5 archive:
            with HDFArchive(self._output_file, 'a') as ar:
                ar[output_group]['iterations'] = iteration_number
                ar[output_group]['chemical_potential'][str(iteration_number)] = self._chemical_potential
                ar[output_group]['dc_imp'][str(iteration_number)] = self._dc_imp
                ar[output_group]['dc_energ'][str(iteration_number)] = self._dc_energ

            # Save the history of Sigma in DCore format
            with h5py.File(self._output_file, 'a') as ar:
                for ish in range(self._n_inequiv_shells):
                    for bname, g in self._sh_quant[ish].Sigma_iw:
                        path = output_group + '/Sigma_iw/ite{}/sh{}/{}'.format(iteration_number, ish, bname)
                        save_giw(ar, path, g)

            # convergence check
            tol = self._params["control"]["converge_tol"]
            if tol > 0:
                print("\nConvergence check")
                x1 = quantities_to_check()
                max_error = numpy.max(numpy.abs(x1 - x0))
                print(" | converge_tol = %.1e" %tol)
                print(" | max_error    = %.1e" %max_error)
                if max_error < tol:
                    print(" | converged --- iteration = %d" %iteration_number)
                    break
                x0 = x1

            sys.stdout.flush()

        self._previous_runs += max_step

    def chemical_potential(self, iteration_number):
        with HDFArchive(self._output_file, 'r') as ar:
            return ar[self._output_group]['chemical_potential'][str(iteration_number)]

    @property
    def n_inequiv_shells(self):
        return self._n_inequiv_shells

    @property
    def inequiv_to_corr(self):
        return self._sk.inequiv_to_corr

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


