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
from itertools import *

from pytriqs.gf.local import *

from dmft_core import DMFTCoreSolver
from matplotlib.gridspec import GridSpec

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
    #

    # Just for convenience
    #output_file = p["model"]["seedname"]+'.out.h5'
    #output_group = 'dmft_out'
    beta = p["system"]["beta"]
    omega_check = p['tool']['omega_check']
    gs = GridSpec(2, 1)

    if fileplot is not None:  # if graph is to be printed in a file
        import matplotlib
        matplotlib.use('Agg')  # do not plot on x11
    from pytriqs.plot.mpl_interface import oplot, plt
    plt.figure(figsize=(8, 10))


    #
    # Load DMFT data
    #
    p['control']['restart'] = True
    solver = DMFTCoreSolver(p["model"]["seedname"], p)
    iteration_number = solver.iteration_number
    nsh = solver.n_inequiv_shells
    spn = solver.spin_block_names
    shell_info = [solver.inequiv_shell_info(ish) for ish in range(nsh)]

    #
    # Chemical potential
    #
    print("  Total number of Iteration: {0}".format(iteration_number))
    print("\n  Iter  Chemical-potential")
    for itr in range(1, iteration_number+1):
        print("  {0} {1}".format(itr, solver.chemical_potential(itr)))

    #
    # Read Sigma and average it
    #
    sigma_ave = []
    nsigma = 0
    num_itr_plot = 7
    itr_sigma = [0]*num_itr_plot
    for itr in range(1, iteration_number+1):
        if itr > iteration_number - num_itr_plot:
            Sigma_iw_sh = solver.Sigma_iw_sh(itr)

            itr_sigma[nsigma] = itr
            sigma_ave.append(GfImFreq(indices=[0], beta=beta, n_points=p["system"]["n_iw"]))
            sigma_ave[nsigma].data[:, 0, 0] = 0.0
            norb_tot = 0
            for ish in range(nsh):
                norb = shell_info[ish]['block_dim']
                for isp in spn:
                    for iorb in range(norb):
                        norb_tot += 1
                        for jorb in range(norb):
                            sigma_ave[nsigma].data[:, 0, 0] += Sigma_iw_sh[ish][isp].data[:, iorb, jorb]
            sigma_ave[nsigma].data[:, 0, 0] /= norb_tot
            nsigma += 1
    #
    # Real part
    #
    plt.subplot(gs[0])
    for itr in range(nsigma):
        oplot(sigma_ave[itr], '-o', mode='R', x_window=(0.0, omega_check), name='Sigma-%s' % itr_sigma[itr])
    plt.legend(loc=0)
    #
    # Imaginary part
    #
    plt.subplot(gs[1])
    for itr in range(nsigma):
        oplot(sigma_ave[itr], '-o', mode='I', x_window=(0.0, omega_check), name='Sigma-%s' % itr_sigma[itr])
    plt.legend(loc=0)

    plt.show()
    if fileplot is not None:
        plt.savefig(fileplot)
    #
    # Output Sigma into a text file
    #
    print("\n Output Local Self Energy : ", p["model"]["seedname"] + "_sigma.dat")
    with open(p["model"]["seedname"] + "_sigma.dat", 'w') as fo:
        print("# Local self energy at imaginary frequency", file=fo)
        #
        # Column information
        #
        print("# [Column] Data", file=fo)
        print("# [1] Frequency", file=fo)
        icol = 1
        for ish in range(nsh):
            norb = shell_info[ish]['block_dim']
            for isp in spn:
                for iorb in range(norb):
                    for jorb in range(norb):
                        icol += 1
                        print("# [%d] Re(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, isp, iorb, jorb), file=fo)
                        icol += 1
                        print("# [%d] Im(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, isp, iorb, jorb), file=fo)
        #
        # Write data
        #
        Sigma_iw_tmp = solver.Sigma_iw_sh(iteration_number)
        omega = [x for x in Sigma_iw_sh[0].mesh]
        for iom in range(len(omega)):
            print("%f " % omega[iom].imag, end="", file=fo)
            for ish in range(nsh):
                norb = shell_info[ish]['block_dim']
                for isp, iorb, jorb in product(spn, range(norb), range(norb)):
                    print("%f %f " % (Sigma_iw_tmp[ish][isp].data[iom, iorb, jorb].real,
                                      Sigma_iw_tmp[ish][isp].data[iom, iorb, jorb].imag), end="", file=fo)
            print("", file=fo)
    #
    # Output Legendre polynomial
    #
    """
    if p["system"]["n_l"] > 0:
        #
        # Output Sigma into a text file
        #
        print("\n Output Local Self Energy : ", p["model"]["seedname"] + "_legendre.dat")
        with open(p["model"]["seedname"] + "_legendre.dat", 'w') as fo:
            print("# Local self energy at imaginary frequency", file=fo)
            #
            # Column information
            #
            print("# [Column] Data", file=fo)
            print("# [1] Order of Legendre polynomials", file=fo)
            icol = 1
            for ish in range(nsh):
                sol[ish].G_l << ar[output_group]['G_l'][str(ish)]
                spn = solver.SK.spin_block_names[solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['SO']]
                norb = solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['dim']
                for isp in spn:
                    for iorb in range(norb):
                        for jorb in range(norb):
                            icol += 1
                            print("# [%d] Re(G_l_{shell=%d, spin=%s, %d, %d})" % (icol, ish, isp, iorb, jorb),
                                  file=fo)
                            icol += 1
                            print("# [%d] Im(G_l_{shell=%d, spin=%s, %d, %d})" % (icol, ish, isp, iorb, jorb),
                                  file=fo)
            #
            # Write data
            #
            for il in range(p["system"]["n_l"]):
                print("%d " % il, end="", file=fo)
                for ish in range(nsh):
                    spn = solver.SK.spin_block_names[solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['SO']]
                    norb = solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['dim']
                    for isp in spn:
                        for iorb in range(norb):
                            for jorb in range(norb):
                                print("%f %f " % (sol[ish].G_l[isp].data[il, iorb, jorb].real,
                                                  sol[ish].G_l[isp].data[il, iorb, jorb].imag), end="", file=fo)
                print("", file=fo)
    """

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
    dcore_check(args.path_input_file, args.output)
