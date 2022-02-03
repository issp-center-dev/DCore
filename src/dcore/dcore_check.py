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

from matplotlib.gridspec import GridSpec
import numpy
import math

from dcore._dispatcher import *
from dcore.dmft_core import DMFTCoreSolver
from dcore.program_options import *

class DMFTCoreCheck(object):

    def __init__(self, ini_file, max_n_iter):
        """
        Parameters
        ----------
        ini_file : string
            Input-file name
        max_n_iter : int
            Max number of iterations to be processed
        """

        if os.path.isfile(ini_file) is False:
            raise Exception("Input file '%s' does not exist." % ini_file)

        print("\n  @ Reading {0} ...".format(ini_file))
        #
        # Construct a parser with default values
        #
        pars = create_parser(['model', 'system', 'tool'])
        #
        # Parse keywords and store
        #
        pars.read(ini_file)
        self.p = pars.as_dict()
        parse_parameters(self.p)

        # Just for convenience
        #output_file = p["model"]["seedname"]+'.out.h5'
        #output_group = 'dmft_out'
        self.beta = self.p["system"]["beta"]

        #
        # Load DMFT data
        #
        self.solver = DMFTCoreSolver(self.p["model"]["seedname"], self.p, read_only=True, restart=True)
        self.n_iter = min(max_n_iter, self.solver.iteration_number)
        self.n_sh = self.solver.n_inequiv_shells
        self.spin_names = self.solver.spin_block_names
        self.shell_info = [self.solver.inequiv_shell_info(ish) for ish in range(self.n_sh)]
        self.n_iw = self.p["system"]["n_iw"]

        print("  Total number of Iteration: {0}".format(self.n_iter))

        # If omega_check is not specified, a fixed number of Matsubara points are taken
        self.omega_check = self.p['tool']['omega_check']
        if self.omega_check == 0:
            nmax = min(30, self.n_iw)
            self.omega_check = (2*nmax+1) * math.pi / self.beta

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

        from ._dispatcher import mpl_interface

        self.plt = mpl_interface.plt
        self.oplot = mpl_interface.oplot


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
                sigma_ave.append(GfImFreq(indices=[0], beta=self.beta, n_points=self.p["system"]["n_iw"], name="$\Sigma_\mathrm{ave}$"))
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
        self.plt.xlim(0, self.omega_check)
        #
        # Imaginary part
        #
        self.plt.subplot(gs[1])
        for itr in range(nsigma):
            self.oplot(sigma_ave[itr], '-o', mode='I', x_window=(0.0, self.omega_check), name='Sigma-%s' % itr_sigma[itr])
        self.plt.legend(loc=0)
        self.plt.xlim(0, self.omega_check)

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

    def _plot_iter(self, basename, fig_ext, data_list, ylabel):
        """
        A core method for plot of convergence

        Parameters
        ----------
        basename
        fig_ext
        data_list: list of dict. The dict must includes
            'y': numpy.ndarray  (mandatory)
            'label': string  ('' if not provided)
            'iorb': int  (0 if not provided)
            'isp': int  (0 if not prvided)
        ylabel: string

        Returns
        -------

        """

        iter = [itr for itr in range(1, self.n_iter + 1)]

        # set linestyles
        linestyles = ['-', ':', '--', '-.', '-', ':']
        markers = ['v', 'x']
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # basic colors
        cmap = self.plt.get_cmap('Set2')
        colors += [cmap(i) for i in range(8)]  # extend colors

        def get_ls(idx):
            return linestyles[idx % len(linestyles)]

        def get_c(idx):
            return colors[idx % len(colors)]

        # plot
        self.__plot_init()
        self.plt.figure(figsize=(8, 10))
        gs = GridSpec(2, 1)

        ax1 = self.plt.subplot(gs[0])
        ax2 = self.plt.subplot(gs[1])

        for data in data_list:
            y = data['y']
            label = data.get('label', '')
            iorb = data.get('iorb', 0)
            isp = data.get('isp', 0)

            ax1.plot(iter, y, label=label, ls=get_ls(isp), marker=markers[isp], color=get_c(iorb))
            diff = abs(numpy.array(y[1:]) - numpy.array(y[:-1]))
            ax2.plot(iter[1:], diff, label=label, ls=get_ls(isp), marker=markers[isp], color=get_c(iorb))

        filename = basename + fig_ext
        ax1.set_xlabel("iterations")
        ax1.set_ylabel(ylabel)
        ax1.set_xlim(left=iter[0])
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=8)

        ax2.set_xlabel("iterations")
        ax2.set_ylabel("diff")
        ax2.set_yscale("log")
        ax2.set_xlim(left=iter[0])
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=8)

        # self.plt.xlim(left=iter[0])
        self.plt.tight_layout()
        self.plt.savefig(filename)
        print(" Output " + filename)

        # save data
        filename = basename + ".dat"
        with open(filename, "w") as f:
            for i, itr in enumerate(iter):
                print(itr, file=f, end="")
                for data in data_list:
                    y = data['y']
                    print("", y[i], file=f, end="")
                print("", file=f)
        print(" Output " + filename)

    def plot_iter_mu(self, basename, fig_ext):
        """
        plot chemical potential as a function of iteration number
        """
        mu = [ self.solver.chemical_potential(itr) for itr in range(1, self.n_iter+1) ]
        data = {'y': mu}
        self._plot_iter(basename, fig_ext, [data, ], "$\mu$")

    def plot_iter_occupation(self, basename, fig_ext):
        """
        plot Occupation number as a function of iteration number
        """

        for ish in range(self.n_sh):
            # Make a graph for each shell
            data_list = []
            norb = self.shell_info[ish]['block_dim']
            for isp, spn in enumerate(self.spin_names):
                for iorb in range(norb):
                    occup = numpy.array([self.solver.density_matrix(itr)[ish][spn][iorb, iorb].real
                                         for itr in range(1, self.n_iter + 1)])

                    data = {
                        'y': occup,
                        'label': "shell=%d, spin=%s, %d" % (ish, spn, iorb),
                        'iorb': iorb,
                        'isp': isp,
                    }
                    data_list.append(data)

            self._plot_iter(basename + '-ish{}'.format(ish), fig_ext, data_list, "Occupation number")

    def plot_iter_spin_moment(self, basename, fig_ext):
        """
        plot spin moment as a function of iteration number
        """

        labels = ["$S_x$", "$S_y$", "$S_z$"]
        for ish in range(self.n_sh):
            # Make a graph for each shell
            data_list = []
            for j in range(3):
                smoment = numpy.array([self.solver.spin_moment(itr)[ish][j]
                                       for itr in range(1, self.n_iter + 1)])

                data = {
                    'y': smoment,
                    'label': "shell=%d, %s" % (ish, labels[j]),
                    'iorb': j,  # change the color depending on the spin direction
                }
                data_list.append(data)

            self._plot_iter(basename + '-ish{}'.format(ish), fig_ext, data_list, "Spin moment")

    def plot_iter_sigma(self, basename, fig_ext):
        """
        plot renormalization factor as a function of iteration number
        """

        w0 = self.n_iw
        for ish in range(self.n_sh):
            # Make a graph for each shell
            data_list = []
            norb = self.shell_info[ish]['block_dim']
            for isp, spn in enumerate(self.spin_names):
                for iorb in range(norb):
                    sigma0 = numpy.array([self.solver.Sigma_iw_sh(itr)[ish][spn].data[w0, iorb, iorb].imag
                                          for itr in range(1, self.n_iter + 1)])
                    z = 1. / (1 - sigma0 / numpy.pi * self.beta)

                    data = {
                        'y': z,
                        'label': "shell=%d, spin=%s, %d" % (ish, spn, iorb),
                        'iorb': iorb,
                        'isp': isp,
                    }
                    data_list.append(data)

            self._plot_iter(basename + '-ish{}'.format(ish), fig_ext, data_list, "Renormalization factor")


def dcore_check(ini_file, prefix, fig_ext, max_n_iter):

    # add a dot to the extension, e.g., 'png' --> '.png'
    ext = fig_ext if fig_ext[0] == '.' else '.' + fig_ext

    # make directory
    dir = os.path.dirname(prefix)
    if not os.path.exists(dir):
        os.makedirs(dir)

    check = DMFTCoreCheck(ini_file, max_n_iter)
    check.print_chemical_potential()
    check.write_sigma_text(basename=prefix+"sigma")
    check.plot_sigma_ave(basename=prefix+"sigma_ave", fig_ext=ext)
    check.plot_iter_mu(basename=prefix+"iter_mu", fig_ext=ext)
    check.plot_iter_sigma(basename=prefix+"iter_sigma", fig_ext=ext)
    check.plot_iter_occupation(basename=prefix+"iter_occup", fig_ext=ext)
    check.plot_iter_spin_moment(basename=prefix+"iter_spin", fig_ext=ext)


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version, print_header

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_check.py',
        description='script for checking the convergence of dcore.',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
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
                        help='prefix for output files (default: check/)'
                        )
    parser.add_argument('--ext',
                        action='store',
                        default='png',
                        type=str,
                        help='file extension of output figures (png, pdf, eps, jpg, etc)'
                        )
    parser.add_argument('--max_n_iter',
                        action='store',
                        default=10000,
                        type=int,
                        help='Max number of iterations to be processed'
                        )
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))
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

    dcore_check(args.path_input_file, args.prefix, args.ext, args.max_n_iter)

    # Finish
    print("\n  Done\n")
