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

import sys, os, copy, re
import __builtin__
import pytriqs.utility.mpi as mpi
from pytriqs.operators.util import *
from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.applications.dft.sumk_dft import *

from typed_parser import *
from program_options import *

def __gettype(name):
    t = getattr(__builtin__, name)
    if isinstance(t, type):
        return t
    raise ValueError(name)


def create_solver_params(dict):
    """
    Parse a dict and create parameters for an impurity solver.
    In dict, keyname should be parameter_name{python_type_name} (e.g. max_time{int}).

    :param dict: a dict object containing keys and values read from *.ini
    :return: a dict object containing parameters
    """
    solver_params = {}
    for k,v in dict.items():
        if k == 'name':
            continue

        r = re.compile('^(.*)\{(.*)\}$')
        try:
            m = r.search(k)
            param_name,param_type_str = m.group(1), m.group(2)
            param_type = __gettype(param_type_str)
        except:
            raise RuntimeError("Unknown type or unrecognized format : " + k)
        solver_params[param_name] = param_type(v)

    return solver_params


class DMFTCoreSolver:
    def __init__(self, seedname, params):
        """
        Initialize solver at each inequivalent correlation shell
        :param seedname:
        :param params:
        """
        self._params = copy.deepcopy(params)
        # Construct a SumKDFT object
        self._SK = SumkDFT(hdf_file=seedname+'.h5', use_dft_blocks=False, h_field=0.0)
        U_file = HDFArchive(seedname+'.h5','r')
        self.Umat = U_file["DCore"]["Umat"]
        self.SO = params['model']['spin_orbit']

        # Construct an impurity solver
        beta = float(params['system']['beta'])
        n_iw = int(params['system']['n_iw']) # Number of Matsubara frequencies
        n_tau = int(params['system']['n_tau']) # Number of tau points
        self._solver_name = params['impurity_solver']['name']

        self._solver_params = create_solver_params(params['impurity_solver'])

        self._h_int = []
        self._S = []
        for ish in range(self._SK.n_inequiv_shells):

            n_orb = self._SK.corr_shells[self._SK.inequiv_to_corr[ish]]['dim']
            if self._SK.SO == 1: n_orb/=2
            spin_names = ["up","down"]
            orb_names = [i for i in range(n_orb)]

            map_operator_structure = {}
            if self.SO:
                for i in range(n_orb):
                    map_operator_structure[('up', i)] = ('ud', i)
                    map_operator_structure[('down', i)] = ('ud', i+n_orb)
            else:
                for i in range(n_orb):
                    map_operator_structure[('up', i)] = ('up', i)
                    map_operator_structure[('down', i)] = ('down', i)

            # Construct Hamiltonian
            self._h_int.append(
                h_int_slater(spin_names, orb_names, U_matrix=self.Umat[self._SK.inequiv_to_corr[ish]], off_diag=True,
                             map_operator_structure=map_operator_structure,
                             H_dump="H" + str(ish) + ".txt", complex=False) # tentative
                )

            # Use GF structure determined by DFT blocks
            gf_struct = self._SK.gf_struct_solver[ish]

            if self._solver_name=="TRIQS/cthyb":
                from pytriqs.applications.impurity_solvers.cthyb import Solver
                self._S.append(Solver(beta=beta, gf_struct=gf_struct, n_iw=n_iw, n_tau=n_tau))
            elif self._solver_name=="TRIQS/hubbard-I":
                from hubbard_solver_matrix import Solver
                self._S.append(Solver(beta=beta, norb=n_orb, use_spin_orbit=self.SO))
            elif self._solver_name=="ALPS/cthyb":
                from pytriqs.applications.impurity_solvers.alps_cthyb import Solver
                self._S.append(Solver(beta=beta, gf_struct=gf_struct, assume_real=True, n_iw=n_iw, n_tau=n_tau))
            else:
                raise RuntimeError("Unknown solver "+self._solver_name)

    # Make read-only getter
    @property
    def Solver(self):
        return self._S

    def solve(self, max_step, output_file, output_group='dmft_output', dry_run=False):
        with_dc = self._params['system']['with_dc']

        fix_mu = self._params['system']['fix_mu']
        if fix_mu:
            mu = self._params['system']['mu']

        sigma_mix = self._params['control']['sigma_mix']  # Mixing factor of Sigma after solution of the AIM
        delta_mix = self._params['control']['delta_mix']  # Mixing factor of Delta as input for the AIM

        prec_mu = self._params['system']['prec_mu']

        previous_runs = 0
        nsh = self._SK.n_inequiv_shells

        # Just for convenience
        SK = self._SK
        S = self._S

        # Set up a HDF file for output
        error = 0
        if mpi.is_master_node():
            try:
                with HDFArchive(output_file, 'a') as f:
                    if output_group in f:
                        if self._params['control']['restart']:
                            ar = f[output_group]
                            if not 'iterations' in ar:
                                raise RuntimeError("Failed to restart the previous simulation!")
    
                            previous_runs = ar['iterations']
                            if ar['iterations'] <= 0:
                                raise RuntimeError("No previous runs to be loaded from " + output_file + "!")
                            print("Loading Sigma_iw... ")
                            for ish in range(nsh):
                                S[ish].Sigma_iw << ar['Sigma_iw'][str(ish)]
                        else:
                            del f[output_group]
                            f.create_group(output_group)
                    else:
                        f.create_group(output_group)
                    f[output_group]['parameters'] = self._params
                    #
                    # Sub group for something
                    #
                    for gname in ['Sigma_iw', 'G0-log', 'G-log', 'Sigma-log', 'G_iw', 'chemical_potential', 'Delta_iw']:
                        if not (gname in f[output_group]): f[output_group].create_group(gname)
            except Exception as e:
                error = 1
                print("Error occurred in IO of a HDF file: " + str(e))
        error = mpi.bcast(error)
        if error != 0:
            return

        previous_runs = mpi.bcast(previous_runs)
        previous_present = previous_runs > 0
        if previous_present:
            if mpi.is_master_node():
                SK.chemical_potential, SK.dc_imp, SK.dc_energ = SK.load(['chemical_potential', 'dc_imp', 'dc_energ'])
                print("Broadcasting Sigma_iw, chemical_potential, dc_imp, dc_energ... ")
            for ish in range(nsh):
                S[ish].Sigma_iw << mpi.bcast(S[ish].Sigma_iw)
                SK.chemical_potential = mpi.bcast(SK.chemical_potential)
                SK.dc_imp = mpi.bcast(SK.dc_imp)
                SK.dc_energ = mpi.bcast(SK.dc_energ)

        for iteration_number in range(previous_runs+1,previous_runs+max_step+1):
            sys.stdout.flush()
            if mpi.is_master_node():
                print("\n########################  Iteration = {0}  ########################\n".format(iteration_number))

            SK.set_Sigma([ S[ish].Sigma_iw for ish in range(nsh) ])                            # set Sigma into the SumK class
            if mpi.is_master_node(): print("  @ ", end='')
            if fix_mu:
                chemical_potential = mu
                chemical_potential = mpi.bcast(chemical_potential)
                SK.set_mu(chemical_potential)
            else:
                SK.calc_mu( precision = prec_mu )  # find the chemical potential for given density
            G_iw_all = SK.extract_G_loc(with_dc=with_dc)
            for ish in range(nsh):
                S[ish].G_iw << G_iw_all[ish]                         # calc the local Green function
                mpi.report("\n    Total charge of Gloc : %.6f"%S[ish].G_iw.total_density())

            # Init the DC term and the real part of Sigma, if no previous runs found:
            if (iteration_number==1 and previous_present==False and with_dc):
                for ish in range(nsh):
                    dm = S[ish].G_iw.density()
                    #SK.calc_dc(dm, orb=ish, U_interact = self._U_int, J_hund = self._J_hund, use_dc_formula = dc_type)
                    self.calc_dc_matrix(dm, orb = ish, Umat=self.Umat[self._SK.inequiv_to_corr[ish]])
                    if self.SO:
                        S[ish].Sigma_iw << SK.dc_imp[self._SK.inequiv_to_corr[ish]]['ud'][0, 0]
                    else:
                        S[ish].Sigma_iw << SK.dc_imp[self._SK.inequiv_to_corr[ish]]['up'][0,0]

            # Calculate new G0_iw to input into the solver:
            for ish in range(nsh):
                if mpi.is_master_node():
                    # We can do a mixing of Delta in order to stabilize the DMFT iterations:
                    S[ish].G0_iw << S[ish].Sigma_iw + inverse(S[ish].G_iw)
                    ar = HDFArchive(output_file, 'a')
                    if (iteration_number>previous_runs+1 or previous_present):
                        mpi.report("\n  @ Mixing input Delta with factor %s"%delta_mix)
                        Delta = (delta_mix * delta(S[ish].G0_iw))\
                                + (1.0-delta_mix) * ar[output_group]['Delta_iw'][str(ish)]
                        S[ish].G0_iw << S[ish].G0_iw + delta(S[ish].G0_iw) - Delta
                    ar[output_group]['Delta_iw'][str(ish)] = delta(S[ish].G0_iw)
                    S[ish].G0_iw << inverse(S[ish].G0_iw)
                    del ar

                S[ish].G0_iw << mpi.bcast(S[ish].G0_iw)

            mpi.report("\n  @ Solve the impurity problem.\n")

            if self._solver_name=="TRIQS/hubbard-I":
                # calculate non-interacting atomic level positions:
                eal = SK.eff_atomic_levels()
                for ish in range(nsh):
                    S[ish].set_atomic_levels( eal = eal[ish] )
                    S[ish].solve(Umat=self.Umat[self._SK.inequiv_to_corr[ish]], verbosity = 0)
            else:
                for ish in range(nsh):
                    S[ish].solve(h_int=self._h_int[ish], **self._solver_params)

            # Solved. Now do post-processing:
            for ish in range(nsh):mpi.report("\n    Total charge of impurity problem : %.6f"%S[ish].G_iw.total_density())

            # Now mix Sigma and G with factor sigma_mix, if wanted:
            if (iteration_number>1 or previous_present):
                if mpi.is_master_node():
                    ar = HDFArchive(output_file,'a')
                    mpi.report("\n  @ Mixing Sigma and G with factor %s"%sigma_mix)
                    for ish in range(nsh):
                        S[ish].Sigma_iw << sigma_mix * S[ish].Sigma_iw \
                                    + (1.0-sigma_mix) * ar[output_group]['Sigma_iw'][str(ish)]
                        S[ish].G_iw << sigma_mix * S[ish].G_iw \
                                + (1.0-sigma_mix) * ar[output_group]['G_iw'][str(ish)]
                    del ar
                for ish in range(nsh):
                    S[ish].G_iw << mpi.bcast(S[ish].G_iw)
                    S[ish].Sigma_iw << mpi.bcast(S[ish].Sigma_iw)

            # Write the final Sigma and G to the hdf5 archive:
            if mpi.is_master_node():
                ar = HDFArchive(output_file, 'a')
                ar[output_group]['iterations'] = iteration_number
                ar[output_group]['chemical_potential'][str(iteration_number)] = SK.chemical_potential
                for ish in range(nsh):
                    ar[output_group]['G_iw'][str(ish)] = S[ish].G_iw
                    ar[output_group]['Sigma_iw'][str(ish)] = S[ish].Sigma_iw
                    #
                    # Save the history of G0, G, Sigma
                    #
                    for gname in ['G0-log', 'G-log', 'Sigma-log']:
                        if not (str(iteration_number) in ar[output_group][gname]):
                            ar[output_group][gname].create_group(str(iteration_number))

                    ar[output_group]['G0-log'][str(iteration_number)][str(ish)] = S[ish].G0_iw
                    ar[output_group]['G-log'][str(iteration_number)][str(ish)] = S[ish].G_iw
                    ar[output_group]['Sigma-log'][str(iteration_number)][str(ish)] = S[ish].Sigma_iw
                del ar

            # Set the new double counting:
            if with_dc:
                for ish in range(nsh):
                    dm = S[ish].G_iw.density() # compute the density matrix of the impurity problem
                    # SK.calc_dc(dm, orb = ish, U_interact = self._U_int, J_hund = self._J_hund, use_dc_formula = dc_type)
                    self.calc_dc_matrix(dm, orb = ish, Umat=self.Umat[SK.inequiv_to_corr[ish]])

            # Save stuff into the user_data group of hdf5 archive in case of rerun:
            SK.save(['chemical_potential','dc_imp','dc_energ'])

    def calc_dc_matrix(self,dens_mat, Umat, orb=0):
        """
        Compute double counting term with U-matrix
        :param dens_mat:gf_struct_solver like
                Density matrix for the specified correlated shell.
        :param Umat: float numpy array [:,:,:,:]
                4-index interaction matrix
        :param orb:int, optional
                Index of an inequivalent shell.
        """
        for icrsh in range(self._SK.n_corr_shells):

            # ish is the index of the inequivalent shell corresponding to icrsh
            ish = self._SK.corr_to_inequiv[icrsh]
            if ish != orb:  continue  # ignore this orbital

            dim = self._SK.corr_shells[icrsh]['dim']

            if self._SK.corr_shells[icrsh]['SO'] == 0:
                spn = self._SK.spin_block_names[self._SK.corr_shells[icrsh]['SO']]
                for sp1 in spn:
                    self._SK.dc_imp[icrsh][sp1] = numpy.zeros((dim,dim), numpy.complex_)
                    for i1 in range(dim):
                        for i2 in range(dim):
                            #
                            # Hartree
                            #
                            for sp2 in spn:
                                self._SK.dc_imp[icrsh][sp1][i1,i2] += numpy.sum(Umat[i1,:,i2,:] * dens_mat[sp2][:,:])
                            #
                            # Exchange
                            #
                            self._SK.dc_imp[icrsh][sp1][i1, i2] += - numpy.sum(Umat[i1, :, :, i2] * dens_mat[sp2][:, :])
            else:
                self._SK.dc_imp[icrsh]["ud"] = numpy.zeros((dim, dim), numpy.complex_)
                dim /= 2
                for s1 in range(2):
                    for i1 in range(dim):
                        for s2 in range(2):
                            for i2 in range(dim):
                                #
                                # Hartree
                                #
                                self._SK.dc_imp[icrsh]["ud"][i1+s1*dim,i2+s1*dim] += numpy.sum(
                                    Umat[i1, :, i2, :] * dens_mat[s2*dim:s2*dim+dim,s2*dim:s2*dim+dim]
                                )
                                #
                                # Exchange
                                #
                                self._SK.dc_imp[icrsh]["ud"][i1 + s1 * dim, i2 + s2 * dim] += numpy.sum(
                                    Umat[i1, :, :, i2] * dens_mat[s2 * dim:s2 * dim + dim, s1 * dim:s1 * dim + dim]
                                )
