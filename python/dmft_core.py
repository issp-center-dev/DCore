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
import re
import time
import __builtin__
import numpy
import copy

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

def is_hermite_conjugate(Sigma_iw):
    """
    Check if Sigma(iw_n) or G(iwn_n) is Hermite conjugate.
    """
    flag = True
    for name, g in Sigma_iw:
        n_points = int(g.data.shape[0]/2)
        for i in range(n_points):
            if not numpy.allclose(g.data[i + n_points, :, :], g.data[n_points - i - 1, :, :].conj().transpose()):
                return False
    return True



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

    def __init__(self, gf_struct, beta, n_iw, n_tau):
        # Self-energy
        self.Sigma_iw = make_block_gf(GfImFreq, gf_struct, beta, n_iw)
        for name, g in self.Sigma_iw:
            g.zero()

        # Lattice Green's function projected onto an impurity
        #self.Gloc_iw = self.Sigma_iw.copy()

        # Impurity Green's function
        #self.Gimp_iw = self.Sigma_iw.copy()


def solve_impurity_model(solver_name, solver_params, mpirun_command, basis_rot, Umat, gf_struct, beta, n_iw, n_tau, Sigma_iw, Gloc_iw, mesh, ish, ite):
    """

    Solve an impurity model

    If mesh is not None, Sigma_w will be computed. Otherwise, None will be returned as Sigma_w.

    """


    Solver = impurity_solvers.solver_classes[solver_name]

    raise_if_mpi_imported()

    sol = Solver(beta, gf_struct, Umat, n_iw, n_tau)

    if not mesh is None and not Solver.is_gf_realomega_available():
        raise RuntimeError("Error: requested real-frequency quantities for an imaginary-time solver.")

    G0_iw = dyson(Sigma_iw=Sigma_iw, G_iw=Gloc_iw)
    sol.set_G0_iw(G0_iw)

    # Compute rotation matrix to the diagonal basis if supported
    rot = compute_diag_basis(G0_iw) if basis_rot else None
    s_params = copy.deepcopy(solver_params)
    s_params['random_seed_offset'] = 1000 * ish
    s_params['work_dir'] = 'work/imp_shell'+str(ish)+"_ite"+str(ite)

    if not mesh is None:
        s_params['calc_Sigma_w'] = True
        s_params['mesh'] = mesh


    # Solve the model
    sol.solve(rot, mpirun_command, s_params)

    # Read & save local quantities
    # Change from DCore v1:
    #      Local impurity Green's function is saved as "Gimp_iw" in DCore v2.
    #      This is intended to distinguish the local impurity Green's function from the one computed by SumkDFT.

    return sol.get_Sigma_iw(), sol.get_Gimp_iw(), sol.get_Sigma_w()



class DMFTCoreSolver(object):
    def __init__(self, seedname, params, output_file='', output_group='dmft_out', read_only=False):
        # Set up
        self._seedname = seedname
        self._params = copy.deepcopy(params)
        self._output_file = seedname+'.out.h5' if output_file is '' else output_file
        self._output_group = output_group
        self._use_spin_orbit = False

        self._beta = float(params['system']['beta'])
        self._n_iw = int(params['system']['n_iw'])  # Number of Matsubara frequencies
        self._n_tau = int(params['system']['n_tau'])  # Number of tau points

        # MPI commands
        if 'num_processes' in params['mpi']:
            self._mpirun_command = params['mpi']['command'].replace('#', str(params['mpi']['num_processes']))

        self._read_only = read_only
        if read_only:
            assert params['control']['restart']

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

        # local quantities at ineuivalent shells
        self._sh_quant = [ShellQuantity(self._gf_struct[ish], self._beta, self._n_iw, self._n_tau) for ish in range(self._n_inequiv_shells)]

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
        if self._params['control']['restart']:
            self._read_output_file__restart()
            assert self._previous_runs >= 1
        else:
            self._prepare_output_file__from_scratch()
            assert self._previous_runs == 0

        self._solver_params = create_solver_params(self._params['impurity_solver'])


        self._sanity_check()


    def _read_output_file__restart(self):
        """
        Read data from & set up an output HDF5 file.
        """

        assert self._params['control']['restart']

        output_file, output_group = self._output_file, self._output_group

        # Set up a HDF file for output
        with HDFArchive(output_file, 'r') as f:
            ar = f[output_group]
            if 'iterations' not in ar:
                raise RuntimeError("Failed to restart the previous simulation!")

            self._previous_runs = ar['iterations']
            if ar['iterations'] <= 0:
                raise RuntimeError("No previous runs to be loaded from " + output_file + "!")

            iterations = ar['iterations']

            print("Loading Sigma_iw... ")
            for ish in range(self._n_inequiv_shells):
                self._sh_quant[ish].Sigma_iw << ar['Sigma_iw'][str(iterations)][str(ish)]

            print("Loading dc_imp and dc_energ... ")
            self._dc_imp = ar['dc_imp'][str(iterations)]
            self._dc_energ = ar['dc_energ'][str(iterations)]
            self._chemical_potential = ar['chemical_potential'][str(iterations)]

    def _prepare_output_file__from_scratch(self):
        """
        Set up an output HDF5 file.
        """

        assert not self._params['control']['restart']

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

        #Gloc_iw_sh, dm_corr_sh = self.calc_Gloc()
        #print("")
        #print("@@@@@@@@@@@@@@@@@@@@@@@@  Double-Counting Correction  @@@@@@@@@@@@@@@@@@@@@@@@")
        #print("")
        #self.calc_dc_imp(dm_corr_sh, set_initial_Sigma_iw=True)

    def _sanity_check(self):
        """
        Sanity checks
        """

        #if not numpy.allclose(self._sk.corr_to_inequiv, self._sk.inequiv_to_corr):
            #raise RuntimeError("corr_to_inequiv must be equal to inequiv_to_corr!")

        raise_if_mpi_imported()

    def _make_sumkdft_params(self):
        return {
            'beta'          : self._params['system']['beta'],
            'prec_mu'       : self._params['system']['prec_mu'],
            'with_dc'       : self._params['system']['with_dc'],
            'Sigma_iw_sh'   : [s.Sigma_iw for s in self._sh_quant],
            'dc_imp'        : self._dc_imp,
            'dc_energ'      : self._dc_energ,
        }

    def calc_Gloc(self):
        """
        Compute the lattice/local Green's function using SumkDFT.

        Return a list of Gloc_iw and density matrices for inequivalent shells
        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'Gloc'
        r = sumkdft.run(self._seedname+'.h5', './work/sumkdft', self._mpirun_command, params)

        if (not self._params['system']['fix_mu']) and (not self._read_only):
            self._chemical_potential = r['mu']

        return r['Gloc_iw_sh'], r['dm_corr_sh']

    def calc_dos(self, Sigma_w_sh, mesh, broadening):
        """

        Compute dos in real frequency.

        :param Sigma_w_sh: list
           List of real-frequency self-energy

        :param broadening: float
           Broadening factor

        :return: tuple
           Results are 'dos', 'dosproj', 'dosproj_orb'.

        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'dos'
        params['mu'] = self._chemical_potential
        params['Sigma_w_sh'] = Sigma_w_sh
        params['mesh'] = mesh
        params['broadening'] = broadening
        r = sumkdft.run(self._seedname+'.h5', './work/sumkdft_dos', self._mpirun_command, params)
        return r['dos'], r['dosproj'], r['dosproj_orb']

    def calc_dos0(self, mesh, broadening):
        """

        Compute dos in real frequency.

        :param broadening: float
           Broadening factor

        :return: tuple
           Results are 'dos0', 'dosproj0', 'dosproj_orb0'.

        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'dos0'
        params['mu'] = self._params['system']['mu']
        params['mesh'] = mesh
        params['broadening'] = broadening
        r = sumkdft.run(self._seedname+'.h5', './work/sumkdft_dos0', self._mpirun_command, params)
        return r['dos0'], r['dosproj0'], r['dosproj_orb0']

    def calc_spaghettis(self, Sigma_w_sh, mesh, broadening):
        """

        Compute A(k, omega)

        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'spaghettis'
        params['mu'] = self._chemical_potential
        params['Sigma_w_sh'] = Sigma_w_sh
        params['mesh'] = mesh
        params['broadening'] = broadening
        r = sumkdft.run(self._seedname+'.h5', './work/sumkdft_spaghettis', self._mpirun_command, params)
        return r['akw']

    def calc_momentum_distribution(self):
        """

        Compute momentum distributions and eigen values of H(k)
        Data are taken from bands_data.

        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'momentum_distribution'
        params['mu'] = self._chemical_potential
        r = sumkdft.run(self._seedname+'.h5', './work/sumkdft_momentum_distribution', self._mpirun_command, params)
        return r['den'], r['ev0']

    def calc_Sigma_w(self, mesh):
        """
        Compute Sigma_w for computing band structure
        For an imaginary-time solver, a list of Nones will be returned.

        :param mesh: (float, float, int)
            real-frequency mesh (min, max, num_points)
        
        :return: list of Sigma_w

        """

        solver_name = self._params['impurity_solver']['name']
        Solver = impurity_solvers.solver_classes[solver_name]
        if Solver.is_gf_realomega_available():
            Gloc_iw_sh, _ = self.calc_Gloc()
            _, _, sigma_w = self.solve_impurity_models(Gloc_iw_sh, -1, mesh)
            return sigma_w
        else:
            return [None] * self.n_inequiv_shells


    def print_density_matrix(self, dm_corr_sh):
        print("\nDensity Matrix")
        for icrsh in range(self._n_corr_shells):
            print("\n  Shell ", icrsh)
            for sp in self._spin_block_names:
                print("\n    Spin ", sp)
                for i1 in range(self._sk.corr_shells[icrsh]['dim']):
                    print("          ", end="")
                    for i2 in range(self._sk.corr_shells[icrsh]['dim']):
                        print("{0:.3f} ".format(dm_corr_sh[icrsh][sp][i1, i2]), end="")
                    print("")


    def solve_impurity_models(self, Gloc_iw_sh, iteration_number, mesh=None):
        """

        Solve impurity models for all inequivalent shells

        :param Gloc_iw_sh:
        :param mesh: (float, float, int)
            (om_min, om_max, n_om)
        :return:
        """

        solver_name = self._params['impurity_solver']['name']
        Sigma_iw_sh = []
        Gimp_iw_sh = []
        Sigma_w_sh = []
        for ish in range(self._n_inequiv_shells):
            print("Solving impurity model for inequivalent shell " + str(ish) + " ...")
            sys.stdout.flush()
            Sigma_iw, Gimp_iw, Sigma_w = solve_impurity_model(solver_name, self._solver_params, self._mpirun_command,
                             self._params["impurity_solver"]["basis_rotation"], self._Umat[ish], self._gf_struct[ish],
                                 self._beta, self._n_iw, self._n_tau,
                                 self._sh_quant[ish].Sigma_iw, Gloc_iw_sh[ish], mesh, ish, iteration_number)
            if not is_hermite_conjugate(Sigma_iw):
                raise RuntimeError("Sigma_iw is not hermite conjugate!")
            if not is_hermite_conjugate(Gimp_iw):
                raise RuntimeError("Gimp_iw is not hermite conjugate!")

            Sigma_iw_sh.append(Sigma_iw)
            if not mesh is None:
                Sigma_w_sh.append(Sigma_w)
            Gimp_iw_sh.append(Gimp_iw)

        # FIXME: DON'T CHANGE THE NUMBER OF RETURNED VALUES. USE DICT INSTEAD.
        if mesh is None:
            return Sigma_iw_sh, Gimp_iw_sh
        else:
            return Sigma_iw_sh, Gimp_iw_sh, Sigma_w_sh

    def calc_dc_imp(self, dm_corr_sh, set_initial_Sigma_iw=True):
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

            dens_mat = dm_corr_sh[self._sk.inequiv_to_corr[ish]]

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
                if set_initial_Sigma_iw:
                    self._sh_quant[ish].Sigma_iw << dc_imp_sh['ud'][0, 0]
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
                if set_initial_Sigma_iw:
                    self._sh_quant[ish].Sigma_iw << dc_imp_sh['up'][0, 0]
                self._dc_imp.append(dc_imp_sh)

            if set_initial_Sigma_iw:
                print("\n      DC Self Energy:")
                for sp1 in self._spin_block_names:
                    print("        Spin {0}".format(sp1))
                    for i1 in range(dim_tot):
                        print("          ", end="")
                        for i2 in range(dim_tot):
                            print("{0:.3f} ".format(self._sh_quant[ish].Sigma_iw[sp1].data[0, i1, i2]), end="")
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

        t0 = time.time()
        for iteration_number in range(self._previous_runs+1, self._previous_runs+max_step+1):
            sys.stdout.flush()
            print("")
            print("#####################################################################")
            print("########################  Iteration = %5d  ########################"%iteration_number)
            print("#####################################################################")
            print("")
            print("")
            print("@@@@@@@@@@@@@@@@@@@@@@@@  Chemical potential and G0_imp  @@@@@@@@@@@@@@@@@@@@@@@@")
            print("")

            # Compute Gloc_iw where the chemical potential is adjusted if needed
            Gloc_iw_sh, dm_corr_sh = self.calc_Gloc()
            self.print_density_matrix(dm_corr_sh)

            for ish in range(self._n_inequiv_shells):
                print("\n    Total charge of Gloc_{shell %d} : %.6f" % (ish, Gloc_iw_sh[ish].total_density()))

            # Compute DC corrections and initial guess to self-energy
            if iteration_number == 1 and not previous_present and with_dc:
                print("")
                print("@@@@@@@@@@@@@@@@@@@@@@@@  Double-Counting Correction  @@@@@@@@@@@@@@@@@@@@@@@@")
                print("")
                self.calc_dc_imp(dm_corr_sh, set_initial_Sigma_iw=True)

            print("\nWall Time : %.1f sec" % (time.time() - t0))

            sys.stdout.flush()
            new_Sigma_iw, new_Gimp_iw = self.solve_impurity_models(Gloc_iw_sh, iteration_number)
            sys.stdout.flush()

            # Solved. Now do post-processing:
            for ish in range(self._n_inequiv_shells):
                print("\nTotal charge of impurity problem : %.6f" % new_Gimp_iw[ish].total_density())

            # Symmetrize over spin components
            if self._params["model"]["time_reversal"]:
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

            # Write data to the hdf5 archive:
            with HDFArchive(self._output_file, 'a') as ar:
                ar[output_group]['iterations'] = iteration_number
                ar[output_group]['chemical_potential'][str(iteration_number)] = self._chemical_potential
                ar[output_group]['dc_imp'][str(iteration_number)] = self._dc_imp
                ar[output_group]['dc_energ'][str(iteration_number)] = self._dc_energ
                for ish in range(self._n_inequiv_shells):
                    #
                    # Save the history of Sigma
                    #
                    if not (str(iteration_number) in ar[output_group]['Sigma_iw']):
                        ar[output_group]['Sigma_iw'].create_group(str(iteration_number))
                    ar[output_group]['Sigma_iw'][str(iteration_number)][str(ish)] = self._sh_quant[ish].Sigma_iw

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
            info['num_orb'] = int(self._dim_sh[ish]/2)
        else:
            info['num_orb'] = self._dim_sh[ish]

        info['block_dim'] = self._dim_sh[ish]

        return info

    def corr_shell_info(self, ish):
        return self.inequiv_shell_info(ish)

    def Sigma_iw_sh(self, iteration_number):
        Sigma_iw_sh = []
        with HDFArchive(self._output_file, 'r') as ar:
            for ish in range(self._n_inequiv_shells):
                Sigma_iw_sh.append(ar[self._output_group]['Sigma_iw'][str(iteration_number)][str(ish)])
        return Sigma_iw_sh

    #@property
    #def sumkdft_compat(self):
        #return self._sk


