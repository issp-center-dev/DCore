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
import os
import sys
import numpy

from dcore.program_options import create_parser, parse_parameters
from dcore.tools import make_empty_dir, launch_mpi_subprocesses
from h5.archive import HDFArchive

def dcore_sparse_bse(filename, np=1, prefix="./"):
    # Construct a parser with default values
    print("\n############  Reading Input File  #################\n")
    print("  Input File Name : ", filename)
    print("")
    pars = create_parser(['model', 'mpi', 'sparse_bse'])

    # Parse keywords and store
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    seedname = p["model"]["seedname"]
    mpirun_command = p['mpi']['command'].replace('#', str(np))

    print("\n############  Reading qsample File  #################\n")
    print("  qsample File Name : ", p['sparse_bse']['qsample'])
    print('')
    with open(p['sparse_bse']['qsample'], 'r') as f:
        num_q = int(f.readline())
        qsample = [], [], []
        for i in range(num_q):
            idx_q, qx, qy, qz = map(int, f.readline().split())
            if idx_q != i:
                raise RuntimeError(f"""File format of {p['sparse_bse']['qsample']} is wrong!""")
            qsample[0].append(qx)
            qsample[1].append(qy)
            qsample[2].append(qz)
        qsample = tuple(map(numpy.stack, qsample))
    
    # Current dir
    cwd_org = os.path.abspath(os.getcwd())

    # Move to work dir
    work_dir = cwd_org + '/work/sparse_bse'
    make_empty_dir(work_dir)
    os.chdir(work_dir)

    # Make input.h5
    with HDFArchive('input.h5', 'w') as h:
        h['qsample'] = qsample
    
    # make output directory
    output_dir = os.path.abspath(f'{cwd_org}/{prefix}')
    print(f'output dir is {output_dir}.')
    if not os.path.exists(output_dir):
        print('Creating output dir...')
        os.makedirs(output_dir)

    commands = [sys.executable, "-m", "dcore.sparse_bse.mpi_main"]
    commands.append('input.h5')
    commands.append(f'{cwd_org}/{seedname}_gk.h5')
    commands.append(f'{output_dir}/{seedname}_chi.h5')
    commands.append(f'--vertex_file={cwd_org}/{seedname}_vertex.h5')
    #if p['sparse_bse']['input_vertex_format'].lower() == 'floc':
        #commands.append(f'--Floc_file={cwd_org}/{seedname}_Floc.h5')
    #elif p['sparse_bse']['input_vertex_format'].lower() == 'g2loc':
        #commands.append(f'--g2loc_file={cwd_org}/{seedname}_g2loc.h5')
    #else:
        #raise ValueError("Invalid input_vertex_format!")
    with open('./output', 'w') as stdout_file:
        launch_mpi_subprocesses(mpirun_command, commands, stdout_file)

    os.chdir(cwd_org)

def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version

    parser = argparse.ArgumentParser(
        prog='dcore_sparse_bse.py',
        description='Post-processing script for dcore (sparse bse).',
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
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))
    parser.add_argument('--prefix',
                        action='store',
                        default='./',
                        type=str,
                        help='prefix for output files (default: post/)'
                        )

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)

    dcore_sparse_bse(args.path_input_file, int(args.np), args.prefix)