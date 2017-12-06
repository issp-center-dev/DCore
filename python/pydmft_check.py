#!/usr/bin/env python
from __future__ import print_function
import os
import argparse
from dmft_core import DMFTCoreSolver, create_parser
from pytriqs.applications.dft.sumk_dft_tools import *
from pytriqs.plot.mpl_interface import oplot, plt

def pydmft_check(filename):
    """
    Main routine for checking convergence

    Parameters
    ----------
    filename : string
        Input-file name
    """
    if mpi.is_master_node(): print("\n  @ Reading {0} ...".format(filename))
    #
    # Construct a parser with default values
    #
    parser = create_parser()
    #
    # Parse keywords and store
    #
    parser.read(filename)
    p = parser.as_dict()
    seedname = p["model"]["seedname"]
    #
    solver = DMFTCoreSolver(p["model"]["seedname"], p)
    #
    beta = float(p['system']['beta'])

    # Just for convenience
    SK = solver._SK
    S = solver._S
    output_file = p["model"]["seedname"]+'.out.h5'
    output_group = 'dmft_out'

    if mpi.is_master_node():
        # Read from HDF file
        ar = HDFArchive(output_file, 'r')
        iteration_number = ar[output_group]['iterations']
        print("  Total number of Iteration: {0}".format(iteration_number))
        print("\n  Iter  Chemical-potential")
        for iter in range(1,iteration_number+1):
            S.Sigma_iw << ar[output_group]['Sigma-%s'%(iter)]
            oplot(S.Sigma_iw["up"], '-o', mode='I', x_window  = (0,20), name = 'Sigma-%s'%(iter))
            print("  {0} {1}".format(iter, ar[output_group]['chemical_potential-%s'%(iter)]))
        del ar

        plt.legend(loc = 4)
        plt.show()
    #
    # Finish
    #
    print("\n  Done\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(\
        prog='pydmft_check.py',\
        description='script for checking the convergence of pydmft.',\
        epilog='end',\
        usage = '$ pydmft_check input',\
        add_help= True)
    parser.add_argument('path_input_file', \
                        action = 'store',\
                        default= None,    \
                        type=str, \
                        help = "input file name."
    )

    args=parser.parse_args()
    if(os.path.isfile(args.path_input_file) is False):
        print("Input file is not exist.")
        sys.exit()
    pydmft_check(args.path_input_file)
