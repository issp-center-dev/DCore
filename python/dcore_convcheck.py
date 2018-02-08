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
from matplotlib.gridspec import GridSpec
import numpy

from program_options import *


def dcore_convcheck(filename, fileplot=None):
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
    output_group = 'dmft_out'
    beta = p["system"]["beta"]
    omega_check = p['tool']['omega_check']
    gs = GridSpec(2, 1)

    if fileplot is not None:  # if graph is to be printed in a file
        import matplotlib
        matplotlib.use('Agg')  # do not plot on x11
    from pytriqs.plot.mpl_interface import oplot, plt
    plt.figure(figsize=(8, 10))
    #
    # Chemical potential
    #
    for ifile in range(p['tool']['nfile']):
        ar = HDFArchive(p["model"]["seedname"]+str(ifile)+'.out.h5', 'r')
        nsh = solver.SK.n_inequiv_shells
        #
        # Read Sigma and average it
        #
        sigma_ave = []
        nsigma = 0
        sigma_ave.append(GfImFreq(indices=[0], beta=beta, n_points=p["system"]["n_iw"]))
        sigma_ave[nsigma].data[:, 0, 0] = 0.0
        norb_tot = 0
        for ish in range(nsh):
            spn = solver.SK.spin_block_names[solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['SO']]
            norb = solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['dim']
            sol[ish].Sigma_iw << ar[output_group]['Sigma-log']["1"][str(ish)]
            for isp in spn:
                for iorb in range(norb):
                    norb_tot += 1
                    for jorb in range(norb):
                        sigma_ave[nsigma].data[:, 0, 0] += sol[ish].Sigma_iw[isp].data[:, iorb, jorb]
            sigma_ave[nsigma].data[:, 0, 0] /= norb_tot
            nsigma += 1
        del ar
        #
        # Real part
        #
        plt.subplot(gs[0])
        for itr in range(nsigma):
            oplot(sigma_ave[itr], '-o', mode='R', x_window=(0.0, omega_check), name='Sigma-%s' % ifile)
        plt.legend(loc=0)
        #
        # Imaginary part
        #
        plt.subplot(gs[1])
        for itr in range(nsigma):
            oplot(sigma_ave[itr], '-o', mode='I', x_window=(0.0, omega_check), name='Sigma-%s' % ifile)
        plt.legend(loc=0)

    plt.show()
    if fileplot is not None:
        plt.savefig(fileplot)
    #
    # Finish
    #
    print("\n  Done\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='dcore_convcheck.py',
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
        dcore_convcheck(args.path_input_file, args.output)
