#!/usr/bin/env python
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
import os
import argparse
from dmft_core import DMFTCoreSolver, create_parser
from pytriqs.applications.dft.sumk_dft_tools import *
from pytriqs.plot.mpl_interface import oplot, plt

def dcore_check(filename):
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
            if iter > iteration_number - 7:
                S.Sigma_iw << ar[output_group]['Sigma-%s'%(iter)]
                oplot(S.Sigma_iw["up"][0,0], '-o', mode='I', x_window  = (p['tool']['omega_min'],p['tool']['omega_max']), name = 'Sigma-%s'%(iter))
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
        prog='dcore_check.py',\
        description='script for checking the convergence of dcore.',\
        epilog='end',\
        usage = '$ dcore_check input',\
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
    dcore_check(args.path_input_file)
