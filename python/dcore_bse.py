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
import sys

from dmft_core import DMFTCoreSolver
from program_options import create_parser

def dcore_bse(filename, np=1):
    """
    Main routine for the BSE post-processing tool

    Parameters
    ----------
    filename : string
        Input-file name
    """

    print("\n############  Reading Input File  #################\n")
    print("  Input File Name : ", filename)

    #
    # Construct a parser with default values
    #
    pars = create_parser()

    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    seedname = p["model"]["seedname"]
    p["mpi"]["num_processes"] = np

    #
    # Load DMFT data
    #
    p['control']['restart'] = True
    solver = DMFTCoreSolver(seedname, p, output_file=seedname+'.out.h5', read_only=True)
    if solver.iteration_number == 0:
        raise RuntimeError("Number of iterations is zero!")
    print("Number of iterations :", solver.iteration_number)

    #
    # Compute data for BSE
    #
    solver.calc_bse()

    #
    # Finish
    #
    print("\n#################  Done  #####################\n")


if __name__ == '__main__':
    import argparse
    from .option_tables import generate_all_description

    parser = argparse.ArgumentParser(
        prog='dcore_bse.py',
        description='Post-processing script for dcore (bse).',
        usage='$ dcore_bse input --np 4',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description(),
        add_help=True)
    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name."
                        )
    parser.add_argument('--np', help='Number of MPI processes', required=True)

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file is not exist.")
        sys.exit(-1)
    dcore_bse(args.path_input_file, int(args.np))
