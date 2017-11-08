from __future__ import print_function
import sys, os
import argparse
from pytriqs.applications.pydmft.typed_parser import TypedParser
from pytriqs.applications.pydmft.dmft_core import DMFTCoreSolver

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

    #
    # Parse keywords and store
    #
    parser.read(filename)
    params = parser.as_dict()

    solver = DMFTCoreSolver(seedname, params)
    solver.plot(output_file=seedname+'.out.h5', output_group='dmft_out')

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
    
