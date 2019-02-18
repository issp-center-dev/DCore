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
import numpy

from program_options import *


class DMFTCoreCheck(object):

    def __init__(self, filename, fileplot=None):
        """
        Main routine for checking convergence

        Parameters
        ----------
        filename : string
            Input-file name

        fileplot : string
            Output file name. File format is determined by the extension (pdf, eps, jpg, etc).
        """

        if os.path.isfile(args.path_input_file) is False:
            raise Exception("Input file '%s' does not exist." % args.path_input_file)

        print("\n  @ Reading {0} ...".format(filename))
        #
        # Construct a parser with default values
        #
        pars = create_parser()
        #
        # Parse keywords and store
        #
        pars.read(filename)
        self.p = pars.as_dict()

        # Just for convenience
        #output_file = p["model"]["seedname"]+'.out.h5'
        #output_group = 'dmft_out'
        self.beta = self.p["system"]["beta"]
        self.omega_check = self.p['tool']['omega_check']

        #
        # Load DMFT data
        #
        self.p['control']['restart'] = True
        self.solver = DMFTCoreSolver(self.p["model"]["seedname"], self.p, read_only=True)
        self.n_iter = self.solver.iteration_number
        self.n_sh = self.solver.n_inequiv_shells
        self.spin_names = self.solver.spin_block_names
        self.shell_info = [self.solver.inequiv_shell_info(ish) for ish in range(self.n_sh)]
        self.n_iw = self.p["system"]["n_iw"]

        print("  Total number of Iteration: {0}".format(self.n_iter))

        # if __plot_init() is called
        self.plot_called = False


    def print_chemical_potential(self):
        """
        print chemical potential
        """

        print("\n  Iter  Chemical-potential")
        for itr in range(1, self.n_iter+1):
            print("  {0} {1}".format(itr, self.solver.chemical_potential(itr)))


    def __plot_init(self):
        if self.plot_called:
            self.plt.clf()
            self.plt.figure(figsize=(8, 6))  # default
            return
        self.plot_called = True

        import matplotlib
        matplotlib.use('Agg')  # do not plot on x11

        from pytriqs.plot.mpl_interface import oplot, plt

        self.plt = plt
        self.oplot = oplot


    def plot_sigma_ave(self, filename):
        """
        plot Sigma(iw) averaged over shell, spin and orbital for last several iterations
        """
        self.__plot_init()
        self.plt.figure(figsize=(8, 10))

        sigma_ave = []
        nsigma = 0
        num_itr_plot = 7
        itr_sigma = [0]*num_itr_plot
        for itr in range(1, self.n_iter+1):
            if itr > self.n_iter - num_itr_plot:
                Sigma_iw_sh = self.solver.Sigma_iw_sh(itr)

                itr_sigma[nsigma] = itr
                sigma_ave.append(GfImFreq(indices=[0], beta=self.beta, n_points=self.p["system"]["n_iw"]))
                sigma_ave[nsigma].data[:, 0, 0] = 0.0
                norb_tot = 0
                for ish in range(self.n_sh):
                    norb = self.shell_info[ish]['block_dim']
                    for isp in self.spin_names:
                        for iorb in range(norb):
                            norb_tot += 1
                            for jorb in range(norb):
                                sigma_ave[nsigma].data[:, 0, 0] += Sigma_iw_sh[ish][isp].data[:, iorb, jorb]
                sigma_ave[nsigma].data[:, 0, 0] /= norb_tot
                nsigma += 1

        gs = GridSpec(2, 1)
        #
        # Real part
        #
        self.plt.subplot(gs[0])
        for itr in range(nsigma):
            self.oplot(sigma_ave[itr], '-o', mode='R', x_window=(0.0, self.omega_check), name='Sigma-%s' % itr_sigma[itr])
        self.plt.legend(loc=0)
        #
        # Imaginary part
        #
        self.plt.subplot(gs[1])
        for itr in range(nsigma):
            self.oplot(sigma_ave[itr], '-o', mode='I', x_window=(0.0, self.omega_check), name='Sigma-%s' % itr_sigma[itr])
        self.plt.legend(loc=0)

        # self.plt.show()
        self.plt.savefig(filename)


    def write_sigma_text(self):
        """
        Output Sigma into a text file
        """

        print("\n Output Local Self Energy : ", self.p["model"]["seedname"] + "_sigma.dat")
        with open(self.p["model"]["seedname"] + "_sigma.dat", 'w') as fo:
            print("# Local self energy at imaginary frequency", file=fo)
            #
            # Column information
            #
            print("# [Column] Data", file=fo)
            print("# [1] Frequency", file=fo)
            icol = 1
            for ish in range(self.n_sh):
                norb = self.shell_info[ish]['block_dim']
                for isp in self.spin_names:
                    for iorb in range(norb):
                        for jorb in range(norb):
                            icol += 1
                            print("# [%d] Re(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, isp, iorb, jorb), file=fo)
                            icol += 1
                            print("# [%d] Im(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, isp, iorb, jorb), file=fo)
            #
            # Write data
            #
            Sigma_iw_tmp = self.solver.Sigma_iw_sh(self.n_iter)
            # omega = [x for x in Sigma_iw_sh[0].mesh]
            omega = [x for x in Sigma_iw_tmp[0].mesh]
            for iom in range(len(omega)):
                print("%f " % omega[iom].imag, end="", file=fo)
                for ish in range(self.n_sh):
                    norb = self.shell_info[ish]['block_dim']
                    for isp, iorb, jorb in product(self.spin_names, range(norb), range(norb)):
                        print("%f %f " % (Sigma_iw_tmp[ish][isp].data[iom, iorb, jorb].real,
                                          Sigma_iw_tmp[ish][isp].data[iom, iorb, jorb].imag), end="", file=fo)
                print("", file=fo)


    def plot_iter_mu(self, filename):
        self.__plot_init()

        iter = [ itr for itr in range(1, self.n_iter+1) ]
        mu = [ self.solver.chemical_potential(itr) for itr in range(1, self.n_iter+1) ]
        self.plt.plot(iter, mu, marker="o")
        self.plt.xlabel("iterations")
        self.plt.ylabel("$\mu$")
        self.plt.savefig(filename)


    def plot_iter_sigma(self, filename):
        self.__plot_init()

        iter = [ itr for itr in range(1, self.n_iter+1) ]

        w0 = self.n_iw
        for ish in range(self.n_sh):
            norb = self.shell_info[ish]['block_dim']
            for isp, iorb, jorb in product(self.spin_names, range(norb), range(norb)):
                sigma0 = numpy.array([ self.solver.Sigma_iw_sh(itr)[ish][isp].data[w0, iorb, jorb].imag
                                       for itr in range(1, self.n_iter+1) ])
                z = 1./(1-sigma0/numpy.pi*self.beta)
                self.plt.plot(iter, z, marker="o", label="(%d,%s,%d,%d)" %(ish,isp,iorb,jorb))

        self.plt.legend()
        self.plt.xlabel("iterations")
        # self.plt.ylabel("$\Sigma(\omega_0)$")
        self.plt.ylabel("Renormalization factor")
        self.plt.savefig(filename)


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

    check = DMFTCoreCheck(args.path_input_file, args.output)
    check.print_chemical_potential()
    check.write_sigma_text()
    check.plot_sigma_ave(args.output)
    check.plot_iter_mu("iter_mu.pdf")
    check.plot_iter_sigma("iter_sigma.pdf")

    # Finish
    print("\n  Done\n")
