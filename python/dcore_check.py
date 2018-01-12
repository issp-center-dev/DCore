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
import argparse
from dmft_core import DMFTCoreSolver
from pytriqs.applications.dft.sumk_dft_tools import *

from program_options import *


def dcore_check(filename, fileplot=None):
    """
    Main routine for checking convergence

    Parameters
    ----------
    filename : string
        Input-file name

    fileplot : string
        Output file name. File format is determined by the extension (pdf, eps, jpg, etc).
    """
    if mpi.is_master_node():
        print("\n  @ Reading {0} ...".format(filename))
    #
    # Construct a parser with default values
    #
    pars = create_parser()
    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    #
    solver = DMFTCoreSolver(p["model"]["seedname"], p)
    #

    # Just for convenience
    sol = solver.S
    output_file = p["model"]["seedname"]+'.out.h5'
    output_group = 'dmft_out'

    if fileplot is not None:  # if graph is to be printed in a file
        import matplotlib
        matplotlib.use('Agg')  # do not plot on x11
    from pytriqs.plot.mpl_interface import oplot, plt

    # Read from HDF file
    ar = HDFArchive(output_file, 'r')
    iteration_number = ar[output_group]['iterations']
    print("  Total number of Iteration: {0}".format(iteration_number))
    print("\n  Iter  Chemical-potential")
    for itr in range(1, iteration_number+1):
        if itr > iteration_number - 7:
            sol[0].Sigma_iw << ar[output_group]['Sigma-log'][str(itr)]['0']
            if solver.SO:
                oplot(sol[0].Sigma_iw["ud"][0, 0], '-o', mode='I',
                      x_window=(p['tool']['omega_min'], p['tool']['omega_max']), name='Sigma-%s' % itr)
            else:
                oplot(sol[0].Sigma_iw["up"][0, 0], '-o', mode='I',
                      x_window=(p['tool']['omega_min'], p['tool']['omega_max']), name='Sigma-%s' % itr)
        print("  {0} {1}".format(itr, ar[output_group]['chemical_potential'][str(itr)]))
    del ar

    plt.legend(loc=4)
    plt.show()
    if fileplot is not None:
        plt.savefig(fileplot)
    #
    # Finish
    #
    print("\n  Done\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='dcore_check.py',
        description='script for checking the convergence of dcore.',
        epilog='end',
        add_help=True)
    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name."
                        )
    parser.add_argument('--output',
                        action='store',
                        default=None,
                        type=str,
                        help='output file name (extension pdf, eps, jpg, etc)'
                        )

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file is not exist.")
        sys.exit(-1)
    if mpi.is_master_node():
        dcore_check(args.path_input_file, args.output)
