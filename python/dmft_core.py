from __future__ import print_function

import sys, os, copy
import pytriqs.utility.mpi as mpi
from pytriqs.operators.util import *
from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.applications.dft.sumk_dft import *


class DMFTCoreSolver:
    def __init__(self, seedname, params):
        self._params = copy.deepcopy(params)

        # Construct a SumKDFT object
        self._SK = SumkDFT(hdf_file=seedname+'.h5', use_dft_blocks=False, h_field=0.0)
        U_file = HDFArchive(seedname+'.h5','r')
        #Umat = U_file["pyDMFT"]["U_matrix"]
        #Upmat = U_file["pyDMFT"]["Up_matrix"]
        self._J_hund = U_file["pyDMFT"]["J_hund"]
        self._U_int = U_file["pyDMFT"]["U_int"]

        # Assume up and down sectors are equivalent. (=paramagnetic)
        self._SK.deg_shells = [[['up','down']]]

        n_orb = self._SK.corr_shells[0]['dim']
        l = self._SK.corr_shells[0]['l']
        spin_names = ["up","down"]
        orb_names = [i for i in range(n_orb)]

        # Construct U matrix for density-density calculations
        Umat, Upmat = U_matrix_kanamori(n_orb=n_orb, U_int=self._U_int, J_hund=self._J_hund)

        # Construct Hamiltonian
        self._h_int = h_int_density(spin_names, orb_names, map_operator_structure=self._SK.sumk_to_solver[0], U=Umat, Uprime=Upmat, H_dump="H.txt")

        # Use GF structure determined by DFT blocks
        gf_struct = self._SK.gf_struct_solver[0]

        # Construct an impurity solver
        beta = float(params['system']['beta'])
        n_iw = int(params['system']['n_iw']) # Number of Matsubara frequencies
        n_tau = int(params['system']['n_tau']) # Number of tau points
        self._name = params['impurity_solver']['name']
        n_l = int(params['impurity_solver']['N_l'])                           # Number of Legendre polynomials
        self._solver_params = {}
        if self._name=="TRIQS/cthyb":
            from pytriqs.applications.impurity_solvers.cthyb import Solver
            self._solver_params["max_time"] = -1
            self._solver_params["random_seed"] = 123 * mpi.rank + 567
            self._solver_params["length_cycle"] = 200
            self._solver_params["n_warmup_cycles"] = 5000
            self._solver_params["n_cycles"] = 50000
            self._S = Solver(beta=beta, gf_struct=gf_struct, n_iw=n_iw, n_tau=n_tau, n_l=n_l)
        elif self._name=="TRIQS/hubbard-I":
            if l==0:
                from hubbard_solver_l0 import Solver
            else:
                from hubbard_solver import Solver
            self._S = Solver(beta=beta, l=l)
        elif self._name=="ALPS/cthyb":
            from pytriqs.applications.impurity_solvers.alps_cthyb import Solver
            self._solver_params["max_time"] = 60     # Max simulation time
            self._solver_params["perform_post_proc"] = True     # Max simulation time
            self._solver_params["verbosity"] = 1     # Max simulation time
            self._S = Solver(beta=beta, gf_struct=gf_struct, assume_real=True, n_l=n_l, n_iw=n_iw, n_tau=n_tau)
        else:
            raise RuntimeError("Unknown solver "+self._name)

    def solve(self, max_step, output_file, output_group='dmft_output', dry_run=False):
        beta = float(self._params['system']['beta'])
        dc_type = int(self._params['system']['dc_type'])                        # DC type: -1 None, 0 FLL, 1 Held, 2 AMF
        if dc_type != -1:
            raise RuntimeError("dc_type != -1 is not supported!")
        fix_mu = self._params['system']['fix_mu']
        if fix_mu:
            mu_init = self._params['system']['mu_init']

        sigma_mix = self._params['control']['sigma_mix']                  # Mixing factor of Sigma after solution of the AIM
        delta_mix = self._params['control']['delta_mix']                  # Mixing factor of Delta as input for the AIM

        prec_mu = 0.0001

        previous_runs = 0
        previous_present = False

        # Just for convenience
        SK = self._SK
        S = self._S

        # Set up a HDF file for output
        if mpi.is_master_node():
            with HDFArchive(output_file, 'a') as f:
                if output_group in f:
                    if self._params['control']['restart']:
                        ar = f[output_group]
                        if 'iterations' in ar:
                            previous_present = True
                            previous_runs = ar['iterations']
                    else:
                        del f[output_group]
                        f.create_group(output_group)
                else:
                    f.create_group(output_group)

        for iteration_number in range(previous_runs+1,previous_runs+max_step+1):
            if mpi.is_master_node():
                print("Iteration = ", iteration_number)

            SK.symm_deg_gf(S.Sigma_iw,orb=0)                        # symmetrise Sigma
            SK.set_Sigma([ S.Sigma_iw ])                            # set Sigma into the SumK class
            if fix_mu:
                chemical_potential = mu_init
                chemical_potential = mpi.bcast(chemical_potential)
                SK.set_mu(chemical_potential)
            else:
                chemical_potential = SK.calc_mu( precision = prec_mu )  # find the chemical potential for given density
            S.G_iw << SK.extract_G_loc()[0]                         # calc the local Green function
            mpi.report("Total charge of Gloc : %.6f"%S.G_iw.total_density())

            # Init the DC term and the real part of Sigma, if no previous runs found:
            if (iteration_number==1 and previous_present==False):
                dm = S.G_iw.density()
                if dc_type >= 0:
                    SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)
                S.Sigma_iw << SK.dc_imp[0]['up'][0,0]

            # Calculate new G0_iw to input into the solver:
            if mpi.is_master_node():
                # We can do a mixing of Delta in order to stabilize the DMFT iterations:
                S.G0_iw << S.Sigma_iw + inverse(S.G_iw)
                ar = HDFArchive(output_file, 'a')
                if (iteration_number>previous_runs+1 or previous_present):
                    mpi.report("Mixing input Delta with factor %s"%delta_mix)
                    Delta = (delta_mix * delta(S.G0_iw)) + (1.0-delta_mix) * ar[output_group]['Delta_iw']
                    S.G0_iw << S.G0_iw + delta(S.G0_iw) - Delta
                ar[output_group]['Delta_iw'] = delta(S.G0_iw)
                S.G0_iw << inverse(S.G0_iw)
                del ar

            #DEBUG
            #print(dir(S.G0_iw['up']))
            #gt = GfImTime(indices = [1, 2, 3], beta = 1.0)
            #gt << InverseFourier(S.G0_iw['up'])
            #print(gt.data)
            S.G0_iw << mpi.bcast(S.G0_iw)

            # Solve the impurity problem:
            if self._name=="TRIQS/hubbard-I":
                S.solve(U_int=self._U_int, J_hund=self._J_hund)
            else:
                S.solve(h_int=self._h_int, **self._solver_params)

            # Solved. Now do post-processing:
            mpi.report("Total charge of impurity problem : %.6f"%S.G_iw.total_density())

            # Now mix Sigma and G with factor sigma_mix, if wanted:
            if (iteration_number>1 or previous_present):
                if mpi.is_master_node():
                    ar = HDFArchive(output_file,'a')
                    mpi.report("Mixing Sigma and G with factor %s"%sigma_mix)
                    S.Sigma_iw << sigma_mix * S.Sigma_iw + (1.0-sigma_mix) * ar[output_group]['Sigma_iw']
                    S.G_iw << sigma_mix * S.G_iw + (1.0-sigma_mix) * ar[output_group]['G_iw']
                    del ar
                S.G_iw << mpi.bcast(S.G_iw)
                S.Sigma_iw << mpi.bcast(S.Sigma_iw)

            # Write the final Sigma and G to the hdf5 archive:
            if mpi.is_master_node():
                ar = HDFArchive(output_file, 'a')
                ar[output_group]['iterations'] = iteration_number + previous_runs
                ar[output_group]['G_iw'] = S.G_iw
                ar[output_group]['Sigma_iw'] = S.Sigma_iw
                ar[output_group]['G0-%s'%(iteration_number + previous_runs)] = S.G0_iw
                ar[output_group]['G-%s'%(iteration_number + previous_runs)] = S.G_iw
                ar[output_group]['Sigma-%s'%(iteration_number + previous_runs)] = S.Sigma_iw
                del ar

            # Set the new double counting:
            dm = S.G_iw.density() # compute the density matrix of the impurity problem
            if dc_type >= 0:
                SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)

            # Save stuff into the user_data group of hdf5 archive in case of rerun:
            SK.save(['chemical_potential','dc_imp','dc_energ'])
