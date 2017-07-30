from __future__ import print_function

import sys, os, copy
import pytriqs.utility.mpi as mpi
from pytriqs.operators.util import *
from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.applications.dft.sumk_dft import *


class DMFTSolver:
    def __init__(self, seedname, params):
        self._params = copy.deepcopy(params)

        # Construct a SumKDFT object
        self._SK = SumkDFT(hdf_file=seedname+'.h5', use_dft_blocks=False, h_field=0.0)
        U_file = HDFArchive(seedname+'.h5','r')
        Umat = U_file["pyDMFT"]["U_matrix"]
        Upmat = U_file["pyDMFT"]["Up_matrix"]

        # Assume up and down sectors are equivalent. (=paramagnetic)
        self._SK.deg_shells = [[['up','down']]]

        n_orb = self._SK.corr_shells[0]['dim']
        l = self._SK.corr_shells[0]['l']
        spin_names = ["up","down"]
        orb_names = [i for i in range(n_orb)]

        # Construct Hamiltonian
        self._h_int = h_int_density(spin_names, orb_names, map_operator_structure=self._SK.sumk_to_solver[0], U=Umat, Uprime=Upmat, H_dump="H.txt")

        # Use GF structure determined by DFT blocks
        gf_struct = self._SK.gf_struct_solver[0]

        # Construct an impurity solver
        beta = params['system']['beta']
        n_iw = params['system']['n_iw'] # Number of Matsubara frequencies
        n_tau = params['system']['n_tau'] # Number of tau points
        solver = params['impurity_solver']['name']
        n_l = params['impurity_solver']['N_l']                            # Number of Legendre polynomials
        self._solver_params = {}
        if solver=="TRIQS/cthyb":
            from pytriqs.applications.impurity_solvers.cthyb import Solver
            self._solver_params["max_time"] = -1
            self._solver_params["random_seed"] = 123 * mpi.rank + 567
            self._solver_params["length_cycle"] = 200
            self._solver_params["n_warmup_cycles"] = 5000
            self._solver_params["n_cycles"] = 50000
            self._S = Solver(beta=beta, gf_struct=gf_struct, n_iw=n_iw, n_tau=n_tau, n_l=n_l)
        elif solver=="ALPS/cthyb":
            from pytriqs.applications.impurity_solvers.alps_cthyb import Solver
            self._solver_params["max_time"] = 60     # Max simulation time
            self._solver_params["perform_post_proc"] = True     # Max simulation time
            self._solver_params["verbosity"] = 1     # Max simulation time
            self._S = Solver(beta=beta, gf_struct=gf_struct, assume_real=True, n_l=n_l, n_iw=n_iw, n_tau=n_tau)
        else:
            raise RuntimeError("Unknown solver "+solver)

    def solve(self, max_step, output_file, output_path):
        beta = self._params['system']['beta']
        dc_type = self._params['system']['.dc_type']                        # DC type: -1 None, 0 FLL, 1 Held, 2 AMF
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

        for iteration_number in range(1,max_step+1):
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
                if (iteration_number>1 or previous_present):
                    mpi.report("Mixing input Delta with factor %s"%delta_mix)
                    Delta = (delta_mix * delta(S.G0_iw)) + (1.0-delta_mix) * ar[output_path]['Delta_iw']
                    S.G0_iw << S.G0_iw + delta(S.G0_iw) - Delta
                ar[output_path]['Delta_iw'] = delta(S.G0_iw)
                S.G0_iw << inverse(S.G0_iw)
                del ar

            S.G0_iw << mpi.bcast(S.G0_iw)

            # Solve the impurity problem:
            S.solve(h_int=self._h_int, **self._solver_params)

            # Solved. Now do post-processing:
            mpi.report("Total charge of impurity problem : %.6f"%S.G_iw.total_density())

            # Now mix Sigma and G with factor sigma_mix, if wanted:
            if (iteration_number>1 or previous_present):
                if mpi.is_master_node():
                    ar = HDFArchive(output_file,'a')
                    mpi.report("Mixing Sigma and G with factor %s"%sigma_mix)
                    S.Sigma_iw << sigma_mix * S.Sigma_iw + (1.0-sigma_mix) * ar[output_path]['Sigma_iw']
                    S.G_iw << sigma_mix * S.G_iw + (1.0-sigma_mix) * ar[output_path]['G_iw']
                    del ar
                S.G_iw << mpi.bcast(S.G_iw)
                S.Sigma_iw << mpi.bcast(S.Sigma_iw)

            # Write the final Sigma and G to the hdf5 archive:
            if mpi.is_master_node():
                ar = HDFArchive(output_file, 'a')
                ar[output_path]['iterations'] = iteration_number + previous_runs
                ar[output_path]['G_iw'] = S.G_iw
                ar[output_path]['Sigma_iw'] = S.Sigma_iw
                ar[output_path]['G0-%s'%(iteration_number + previous_runs)] = S.G0_iw
                ar[output_path]['G-%s'%(iteration_number + previous_runs)] = S.G_iw
                ar[output_path]['Sigma-%s'%(iteration_number + previous_runs)] = S.Sigma_iw
                del ar

            # Set the new double counting:
            dm = S.G_iw.density() # compute the density matrix of the impurity problem
            if dc_type >= 0:
                SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)

            # Save stuff into the user_data group of hdf5 archive in case of rerun:
            SK.save(['chemical_potential','dc_imp','dc_energ'])
