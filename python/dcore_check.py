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
    ar = HDFArchive(output_file, 'r')
    iteration_number = ar[output_group]['iterations']
    perform_tail_fit = p["system"]["perform_tail_fit"] and \
                       not ar[output_group]['parameters']["system"]["perform_tail_fit"]
    nsh = solver.SK.n_inequiv_shells
    print("  Total number of Iteration: {0}".format(iteration_number))
    print("\n  Iter  Chemical-potential")
    for itr in range(1, iteration_number+1):
        print("  {0} {1}".format(itr, ar[output_group]['chemical_potential'][str(itr)]))
    #
    # Read Sigma and average it
    #
    sigma_ave = []
    sigma_fit = []
    nsigma = 0
    itr_sigma = [0]*7
    for itr in range(1, iteration_number+1):
        if itr > iteration_number - 7:
            itr_sigma[nsigma] = itr
            sigma_ave.append(GfImFreq(indices=[0], beta=beta, n_points=p["system"]["n_iw"]))
            sigma_fit.append(GfImFreq(indices=[0], beta=beta, n_points=p["system"]["n_iw"]))
            sigma_ave[nsigma].data[:, 0, 0] = 0.0
            sigma_fit[nsigma].data[:, 0, 0] = 0.0
            norb_tot = 0
            for ish in range(nsh):
                spn = solver.SK.spin_block_names[solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['SO']]
                norb = solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['dim']
                sol[ish].Sigma_iw << ar[output_group]['Sigma-log'][str(itr)][str(ish)]
                sigma_iw_fit = sol[ish].Sigma_iw.copy()
                if perform_tail_fit:
                    tail_fit(sigma_iw_fit, fit_max_moment=p["system"]["fit_max_moment"],
                             fit_min_w=p["system"]["fit_min_w"], fit_max_w=p["system"]["fit_max_w"])
                for isp in spn:
                    for iorb in range(norb):
                        norb_tot += 1
                        for jorb in range(norb):
                            sigma_ave[nsigma].data[:, 0, 0] += sol[ish].Sigma_iw[isp].data[:, iorb, jorb]
                            sigma_fit[nsigma].data[:, 0, 0] += sigma_iw_fit[isp].data[:, iorb, jorb]
            sigma_ave[nsigma].data[:, 0, 0] /= norb_tot
            sigma_fit[nsigma].data[:, 0, 0] /= norb_tot
            nsigma += 1
    del ar
    #
    # Real part
    #
    plt.subplot(gs[0])
    for itr in range(nsigma):
        oplot(sigma_ave[itr], '-o', mode='R', x_window=(0.0, omega_check), name='Sigma-%s' % itr_sigma[itr])
        if perform_tail_fit:
            oplot(sigma_fit[itr], '-o', mode='R', x_window=(0.0, omega_check), name='S_fit-%s' % itr_sigma[itr])
    plt.legend(loc=0)
    #
    # Imaginary part
    #
    plt.subplot(gs[1])
    for itr in range(nsigma):
        oplot(sigma_ave[itr], '-o', mode='I', x_window=(0.0, omega_check), name='Sigma-%s' % itr_sigma[itr])
        if perform_tail_fit:
            oplot(sigma_fit[itr], '-o', mode='I', x_window=(0.0, omega_check), name='S_fit-%s' % itr_sigma[itr])
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
            spn = solver.SK.spin_block_names[solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['SO']]
            norb = solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['dim']
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
        omega = [x for x in sol[0].Sigma_iw.mesh]
        for iom in range(len(omega)):
            print("%f " % omega[iom].imag, end="", file=fo)
            for ish in range(nsh):
                spn = solver.SK.spin_block_names[solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['SO']]
                norb = solver.SK.corr_shells[solver.SK.inequiv_to_corr[ish]]['dim']
                for isp in spn:
                    for iorb in range(norb):
                        for jorb in range(norb):
                            print("%f %f " % (sol[ish].Sigma_iw[isp].data[iom, iorb, jorb].real,
                                              sol[ish].Sigma_iw[isp].data[iom, iorb, jorb].imag), end="", file=fo)
            print("", file=fo)
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
