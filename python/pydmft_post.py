from __future__ import print_function
import sys, os, copy
import argparse
import pytriqs.utility.mpi as mpi
from pytriqs.archive.hdf_archive import HDFArchive
from pytriqs.applications.pydmft.typed_parser import TypedParser
from pytriqs.applications.pydmft.dmft_core import DMFTCoreSolver
from pytriqs.applications.dft.sumk_dft import *

def pydmft_post(filename, seedname):

    #
    # Set Default value
    #
    parser = TypedParser()
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature")
    parser.add_option("system", "n_iw", int, 2048, "Number of Matsubara frequencies")
    parser.add_option("system", "n_tau", int, 10000, "Number of imaginary-time points")
    parser.add_option("system", "dc_type", int, -1, "Type of double-counting correction")
    parser.add_option("system", "fix_mu", bool, False, "Whether or not to use a fixed chemical potential")
    
    parser.add_option("impurity_solver", "N_l", int, 50, "Number of Legendre polynomials")
    parser.add_option("impurity_solver", "name", str, 'TRIQS/cthyb', "Name of impurity solver")
    
    parser.add_option("control", "max_step", int, 100, "Max number of SCF steps")
    parser.add_option("control", "sigma_mix", float, 0.5, "Mixing parameter for self-energy")
    parser.add_option("control", "delta_mix", float, 0.5, "Mixing parameter for hybridization function")
    parser.add_option("control", "restart", bool, False, "Whether or not restart from a previous calculation")

    parser.add_option("dos", "flag", bool, False, "Whether or not calculate DOS")
    parser.add_option("band", "flag", bool, False, "Whether or not calculate band")

    # ToDo set default values
    parser.add_option("tool", "omega_min", float, 0, "Minimum value of real frequency")
    parser.add_option("tool", "omega_max", float, 0, "Max value of real frequency")
    parser.add_option("tool", "Nomega", int, 0, "Number of real frequencies")
    parser.add_option("tool", "broadening", float, 0, "An additional Lorentzian broadening")



    #
    # Parse keywords and store
    #
    parser.read(filename)
    params = parser.as_dict()

    dct=DMFTCoreTools(seedname, params)
    dct.post()


    #solver.plot(output_file=seedname+'.out.h5', output_group='dmft_out')

class DMFTCoreTools:
    def __init__(self, seedname, params):
        self._params = copy.deepcopy(params)
        # Construct a SumKDFT object
        self._flag_dos=bool(params['dos']['flag'])
        self._flag_band=bool(params['band']['flag'])
        self._omega_min=float(params['tool']['omega_min'])
        self._omega_max=float(params['tool']['omega_max'])
        self._Nomega=int(params['tool']['Nomega'])
        self._broadning=float(params['tool']['bradning'])
        self._seedname=seedname
        self._SKT = SumkDFTTools(hdf_file=self._seedname + '.h5', use_dft_blocks=False)
        solver = DMFTCoreSolver(seedname, params)
        self._S=solver._S

    def post(self):
        flag_dos=self._flag_dos
        flag_band=self._flag_band
        SKT = self._SKT
        S=self._S
        if flag_dos or flag_band:
            if mpi.is_master_node():
                print("~ Real frequency")

            # set necessary quantities
            if mpi.is_master_node():
                SKT.chemical_potential, SKT.dc_imp, SKT.dc_energ = SKT.load(['chemical_potential','dc_imp','dc_energ'])
            SKT.chemical_potential = mpi.bcast(SKT.chemical_potential)
            SKT.dc_imp = mpi.bcast(SKT.dc_imp)
            SKT.dc_energ = mpi.bcast(SKT.dc_energ)

            if S._name == 'TRIQS/hubbard-I':
                # set atomic levels:
                eal = SKT.eff_atomic_levels()[0]
                S.set_atomic_levels( eal = eal )
                # Run the solver to get GF and self-energy on the real axis
                S.GF_realomega(ommin=self._omega_min, ommax = self._omega_max, N_om=self._Nomega, U_int=S._U_int, J_hund=S._J_hund)
            elif S._name == "TRIQS/cthyb" or S._name == "ALPS/cthyb":
                print("not supported yet")
                return
            else:
                raise RuntimeError("Unknown solver " + S._name)

            SKT.set_Sigma([S.Sigma])

        if flag_dos:
            if mpi.is_master_node():
                # ToDo Time stamp
                print ("~ calc DOS")
            SKT.dos_parproj_basis(broadening=self._broadening)

        if flag_band:
            if mpi.is_master_node():
                # ToDo Time stamp
                print ("~ calc band")
            SKT.spaghettis(broadening=self._broadening,plot_range=None,ishell=None,save_to_file='Akw_')




if __name__ == '__main__':

    parser = argparse.ArgumentParser(\
                                 prog='pydmft_post.py',\
                                 description='.',\
                                 epilog='end',\
                                 usage = '$ pydmft_post input.ini seedname',\
                                 add_help= True)

    parser.add_argument('path_input_file', \
                    action = 'store',\
                    default= None,    \
                    type=str, \
                    help = "input file name.")

    parser.add_argument('seedname', \
                    action = 'store',\
                    default= None,    \
                    type=str, \
                    help = "seed name.")

    args=parser.parse_args()
    if(os.path.isfile(args.path_input_file) is False):
        print("Input file is not exist.")
        sys.exit()

    pydmft_post(args.path_input_file, args.seedname)


