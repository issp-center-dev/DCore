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

from itertools import *

from pytriqs.gf import *

from dmft_core import DMFTCoreSolver
from matplotlib.gridspec import GridSpec
import numpy

from program_options import *


class DMFTCoreCheck(object):

    def __init__(self, ini_file):
        """
        Parameters
        ----------
        ini_file : string
            Input-file name
        """

        if os.path.isfile(ini_file) is False:
            raise Exception("Input file '%s' does not exist." % ini_file)

        print("\n  @ Reading {0} ...".format(ini_file))
        #
        # Construct a parser with default values
        #
        pars = create_parser()
        #
        # Parse keywords and store
        #
        pars.read(ini_file)
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


    def plot_sigma_ave(self, basename, fig_ext):
        """
        plot Sigma(iw) averaged over shell, spin and orbital for last several iterations
        """
        self.__plot_init()

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

        self.plt.figure(figsize=(8, 10))
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

        filename = basename + fig_ext
        self.plt.savefig(filename)
        print(" Output " + filename)


    def write_sigma_text(self, basename):
        """
        Output Sigma into a text file
        """

        from .tools import save_Sigma_iw_sh_txt

        filename = basename + '.dat'
        save_Sigma_iw_sh_txt(basename + '.dat', self.solver.Sigma_iw_sh(self.n_iter), self.spin_names)

        print(" Output " + filename)


    def plot_iter_mu(self, basename, fig_ext):
        """
        plot chemical potential as a function of iteration number
        """
        self.__plot_init()

        iter = [ itr for itr in range(1, self.n_iter+1) ]
        mu = [ self.solver.chemical_potential(itr) for itr in range(1, self.n_iter+1) ]

        # save data
        filename = basename + ".dat"
        numpy.savetxt(filename, zip(iter, mu), fmt="%d %.10e")
        print(" Output " + filename)

        # plot
        filename = basename + fig_ext
        self.plt.figure(figsize=(8, 10))
        gs = GridSpec(2, 1)

        self.plt.subplot(gs[0])
        self.plt.plot(iter, mu, marker="o")
        self.plt.xlabel("iterations")
        self.plt.ylabel("$\mu$")

        # diff
        self.plt.subplot(gs[1])
        self.plt.plot(iter[1:], abs(numpy.array(mu[1:])-numpy.array(mu[:-1])), marker="o")
        self.plt.xlabel("iterations")
        self.plt.ylabel("diff")
        self.plt.yscale("log")

        self.plt.savefig(filename)
        print(" Output " + filename)

    def plot_iter_sigma(self, basename, fig_ext):
        """
        plot renormalization factor as a function of iteration number
        """
        self.__plot_init()

        iter = [ itr for itr in range(1, self.n_iter+1) ]
        label_all = []
        z_all = []

        w0 = self.n_iw
        for ish in range(self.n_sh):
            norb = self.shell_info[ish]['block_dim']
            for isp, iorb, jorb in product(self.spin_names, range(norb), range(norb)):
                sigma0 = numpy.array([ self.solver.Sigma_iw_sh(itr)[ish][isp].data[w0, iorb, jorb].imag
                                       for itr in range(1, self.n_iter+1) ])
                z = 1./(1-sigma0/numpy.pi*self.beta)
                z_all.append(z)
                label_all.append("shell=%d, spin=%s, %d, %d" %(ish,isp,iorb,jorb))

        # save data
        filename = basename + ".dat"
        with open(filename, "w") as f:
            for i, itr in enumerate(iter):
                print(itr, file=f, end="")
                for z in z_all:
                    print("", z[i], file=f, end="")
                print("", file=f)
        print(" Output " + filename)

        # plot
        filename = basename + fig_ext
        self.plt.figure(figsize=(8, 10))
        gs = GridSpec(2, 1)

        self.plt.subplot(gs[0])
        for z, label in zip(z_all, label_all):
            self.plt.plot(iter, z, label=label, marker="o")
        self.plt.xlabel("iterations")
        self.plt.ylabel("Renormalization factor")
        self.plt.legend()

        # diff
        self.plt.subplot(gs[1])
        for z, label in zip(z_all, label_all):
            self.plt.plot(iter[1:], abs(numpy.array(z[1:])-numpy.array(z[:-1])), label=label, marker="o")
        self.plt.xlabel("iterations")
        self.plt.ylabel("diff")
        self.plt.yscale("log")
        self.plt.legend()

        self.plt.savefig(filename)
        print(" Output " + filename)


def dcore_check(ini_file, prefix, fig_ext):

    # add a dot to the extension, e.g., 'png' --> '.png'
    ext = fig_ext if fig_ext[0] == '.' else '.' + fig_ext

    # make directory
    dir = os.path.dirname(prefix)
    if not os.path.exists(dir):
        os.makedirs(dir)

    check = DMFTCoreCheck(ini_file)
    check.print_chemical_potential()
    check.write_sigma_text(basename=prefix+"sigma")
    check.plot_sigma_ave(basename=prefix+"sigma_ave", fig_ext=ext)
    check.plot_iter_mu(basename=prefix+"iter_mu", fig_ext=ext)
    check.plot_iter_sigma(basename=prefix+"iter_sigma", fig_ext=ext)


if __name__ == '__main__':
    from .option_tables import generate_all_description
    import argparse

    parser = argparse.ArgumentParser(
        prog='dcore_check.py',
        description='script for checking the convergence of dcore.',
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
    parser.add_argument('--prefix',
                        action='store',
                        default='check/',
                        type=str,
                        help='prefix for output files'
                        )
    parser.add_argument('--ext',
                        action='store',
                        default='png',
                        type=str,
                        help='file extension of output figures (png, pdf, eps, jpg, etc)'
                        )
    # for backward compatibility
    # parser.add_argument('--output',
    #                     action='store',
    #                     default=None,
    #                     type=str,
    #                     help='not used (retained for backward compatibility)'
    #                     )

    # if args.prefix is not None:
    #     warn("--output option is not used")

    args = parser.parse_args()

    dcore_check(args.path_input_file, args.prefix, args.ext)

    # Finish
    print("\n  Done\n")
