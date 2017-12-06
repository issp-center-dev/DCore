from __future__ import print_function

import sys, os, copy, re
import __builtin__
import pytriqs.utility.mpi as mpi
from pytriqs.operators.util import *
from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.applications.dft.sumk_dft import *

from typed_parser import *


def create_parser():
    parser = TypedParser()

    parser.add_option("model", "t", float, 1.0, "Nearest neighbor hopping")
    parser.add_option("model", "t'", float, 0.0, "Second nearest neighbor hopping")
    parser.add_option("model", "lattice", str, "chain", "Lattice name")
    parser.add_option("model", "orbital_model", str, "single", "some help message")
    parser.add_option("model", "ncor", int, 1, "Number of correlation shell")
    parser.add_option("model", "cshell", str, "[]", "Angular momentum and number of states of each shell")
    parser.add_option("model", "seedname", str, "pydmft", "some help message")
    parser.add_option("model", "U", float, 0.0, "Coulomb")
    parser.add_option("model", "J", float, 0.0, "Hund")
    parser.add_option("model", "nelec", float, 1.0, "Number of electrons")
    parser.add_option("model", "bvec", str, "[(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)]", "Reciprocal lattice vectors")

    parser.add_option("system", "beta", float, 1.0, "Inverse temperature")
    parser.add_option("system", "n_iw", int, 2048, "Number of Matsubara frequencies")
    parser.add_option("system", "n_tau", int, 10000, "Number of imaginary-time points")
    parser.add_option("system", "dc_type", int, -1, "Type of double-counting correction. Set -1 to disable DC correction. Other availale options are according to TRIQS DFTTools.")
    parser.add_option("system", "fix_mu", bool, False, "Whether or not to use a fixed chemical potential")
    parser.add_option("system", "nk", int, 8, "Number of k along each line")
    parser.add_option("system", "nk0", int, 0, "Number of k along b_0")
    parser.add_option("system", "nk1", int, 0, "Number of k along b_1")
    parser.add_option("system", "nk2", int, 0, "Number of k along b_2")

    parser.add_option("impurity_solver", "name", str, 'TRIQS/cthyb', "Name of impurity solver. Available options are TRIQS/cthyb, TRIQS/hubbard-I, ALPS/cthyb.")

    # Allow undefined options
    parser.allow_undefined_options("impurity_solver")

    parser.add_option("control", "max_step", int, 100, "Max number of SCF steps")
    parser.add_option("control", "sigma_mix", float, 0.5, "Mixing parameter for self-energy")
    parser.add_option("control", "delta_mix", float, 0.5, "Mixing parameter for hybridization function")
    parser.add_option("control", "restart", bool, False, "Whether or not restart from a previous calculation")

    parser.add_option("tool", "omega_min", float, -1, "Minimum value of real frequency")
    parser.add_option("tool", "omega_max", float, 1, "Max value of real frequency")
    parser.add_option("tool", "Nomega", int, 100, "Number of real frequencies")
    parser.add_option("tool", "broadening", float, 0.1, "An additional Lorentzian broadening")
    parser.add_option("tool", "eta", float, 0.01, "Imaginary frequency shift")
    parser.add_option("tool", "nnode", int, 2, "Number of k-node for band path")
    parser.add_option("tool", "knode", str, "[(G,0.0,0.0,0.0),(X,1.0,0.0,0.0)]", "k-node for band path")
    parser.add_option("tool", "nk_line", int, 8, "Number of k along each line")

    return parser


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
        self._params = copy.deepcopy(params)
        # Construct a SumKDFT object
        self._SK = SumkDFT(hdf_file=seedname+'.h5', use_dft_blocks=False, h_field=0.0)
        U_file = HDFArchive(seedname+'.h5','r')
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
        self._h_int = h_int_kanamori(spin_names, orb_names, map_operator_structure=self._SK.sumk_to_solver[0], U=Umat, Uprime=Upmat, J_hund=self._J_hund, H_dump="H.txt")

        # Use GF structure determined by DFT blocks
        gf_struct = self._SK.gf_struct_solver[0]

        # Construct an impurity solver
        beta = float(params['system']['beta'])
        n_iw = int(params['system']['n_iw']) # Number of Matsubara frequencies
        n_tau = int(params['system']['n_tau']) # Number of tau points
        self._solver_name = params['impurity_solver']['name']

        self._solver_params = create_solver_params(params['impurity_solver'])

        if self._solver_name=="TRIQS/cthyb":
            from pytriqs.applications.impurity_solvers.cthyb import Solver
            self._S = Solver(beta=beta, gf_struct=gf_struct, n_iw=n_iw, n_tau=n_tau)
        elif self._solver_name=="TRIQS/hubbard-I":
            if l==0:
                from hubbard_solver_l0 import Solver
            else:
                from hubbard_solver import Solver
            self._S = Solver(beta=beta, l=l)
        elif self._solver_name=="ALPS/cthyb":
            from pytriqs.applications.impurity_solvers.alps_cthyb import Solver
            self._S = Solver(beta=beta, gf_struct=gf_struct, assume_real=True, n_iw=n_iw, n_tau=n_tau)
        else:
            raise RuntimeError("Unknown solver "+self._solver_name)

    # Make read-only getter
    @property
    def Solver(self):
        return self._S

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

        # Just for convenience
        SK = self._SK
        S = self._S

        # Set up a HDF file for output
        if mpi.is_master_node():
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
                        S.Sigma_iw << ar['Sigma_iw']
                    else:
                        del f[output_group]
                        f.create_group(output_group)
                else:
                    f.create_group(output_group)

        previous_runs = mpi.bcast(previous_runs)
        previous_present = previous_runs > 0
        if previous_present:
            if mpi.is_master_node():
                print("Broadcasting Sigma_iw... ")
            S.Sigma_iw << mpi.bcast(S.Sigma_iw)

        if mpi.is_master_node():
            ar = HDFArchive(output_file, 'a')
            ar[output_group]['parameters'] = self._params

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

            S.G0_iw << mpi.bcast(S.G0_iw)

            # Solve the impurity problem:
            if self._solver_name=="TRIQS/hubbard-I":
                # calculate non-interacting atomic level positions:
                eal = SK.eff_atomic_levels()[0]
                S.set_atomic_levels( eal = eal )
                S.solve(U_int=self._U_int, J_hund=self._J_hund, verbosity = 1)
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
                ar[output_group]['iterations'] = iteration_number 
                ar[output_group]['G_iw'] = S.G_iw
                ar[output_group]['Sigma_iw'] = S.Sigma_iw
                ar[output_group]['G0-%s'%(iteration_number)] = S.G0_iw
                ar[output_group]['G-%s'%(iteration_number)] = S.G_iw
                ar[output_group]['Sigma-%s'%(iteration_number)] = S.Sigma_iw
                ar[output_group]['chemical_potential-%s'%(iteration_number)] = SK.chemical_potential
                del ar

            # Set the new double counting:
            dm = S.G_iw.density() # compute the density matrix of the impurity problem
            if dc_type >= 0:
                SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)

            # Save stuff into the user_data group of hdf5 archive in case of rerun:
            SK.save(['chemical_potential','dc_imp','dc_energ'])
