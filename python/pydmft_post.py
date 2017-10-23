from __future__ import print_function
import sys, os
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

    #
    # If input file is not specified ... 
    #
    if len(sys.argv) != 3:
        print("Usage:")
        print("$ pydmft_post input seedname")
        sys.exit()
    pydmft_post(sys.argv[1], sys.argv[2])
