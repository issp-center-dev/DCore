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

import sys
from dcore.dmft_core import DMFTCoreSolver
from dcore.program_options import *


def dcore(filename, np=1):
    """
    Main routine of DCore

    Parameters
    ----------
    filename : string
        Input-file name
    """
    # Set Default value
    pars = create_parser(['model', 'system', 'impurity_solver', 'control', 'mpi'])
    #
    # Parse keywords and store
    #
    pars.read(filename)
    params = pars.as_dict()
    parse_parameters(params)

    params["mpi"]["num_processes"] = np

    solver = DMFTCoreSolver(params["model"]["seedname"], params, restart=params['control']['restart'])

    solver.do_steps(max_step=params["control"]["max_step"])

    print("\n########################  Done  ########################\n")


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version, print_header

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore.py',
        description='.',
        usage='$ dcore input.ini --np 4',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description()
    )

    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name.")
    parser.add_argument('--np', default=1, help='Number of MPI processes', required=True)
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)

    dcore(args.path_input_file, int(args.np))
