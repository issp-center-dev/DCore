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

from dcore.program_options import create_parser, parse_parameters
from dcore.tools import launch_mpi_subprocesses
import os

def dcore_mpicheck(filename, np=1):
    """
    Main routine for the dcore_mpicheck

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
    pars = create_parser(['mpi'])
    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    p["mpi"]["num_processes"] = np
    mpirun_command = p['mpi']['command'].replace('#', str(p['mpi']['num_processes']))


    #
    # Summary of input parameters
    #
    print("\n  @ Parameter summary")
    print("\n    [mpi] block")
    for k, v in list(p["mpi"].items()):
        print("      {0} = {1}".format(k, v))

    print("")
    if 'DCORE_MPIRUN_COMMAND' in os.environ:
        print("Environment variable \"DCORE_MPIRUN_COMMAND\" is set to ", os['DCORE_MPIRUN_COMMAND'])
        print("This is used as the default MPI commpand, which can be overitten by the input parameter [mpi][command].")
    else:
        print("Environment variable \"DCORE_MPIRUN_COMMAND\" is not set.")

    print("Actual MPI commmand: ", mpirun_command)
    output_fname = "output_mpicheck.txt"
    print("Output file: ", output_fname)

    other_commands = ["echo", "Hello"]
    print("Calling: ", " ".join([mpirun_command] + other_commands))
    with open(output_fname, 'w') as output_f:
        launch_mpi_subprocesses(mpirun_command, other_commands, output_f)
    print("Done")
    print("Reading output file: ", output_fname)
    with open(output_fname, 'r') as output_f:
        for line in output_f.readlines():
            print(line, end='')
    print("\n#################  Done  #####################\n")


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version, print_header
    import os, sys

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_mpicheck.py',
        description='Checker for MPI environment and input parameterrs',
        usage='$ dcore_mpicheck input --np 4',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description()
    )
    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name."
                        )
    parser.add_argument('--np', help='Number of MPI processes', required=True)
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))
    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)
    dcore_mpicheck(args.path_input_file, int(args.np))
