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
import os
import sys
import numpy
import copy

from pytriqs.archive import HDFArchive
from .pytriqs_gf_compat import *
from pytriqs.operators import *

from dmft_core import DMFTCoreSolver
from sumkdft import SumkDFTCompat
from program_options import create_parser, parse_parameters

from .tools import launch_mpi_subprocesses
import impurity_solvers
from . import sumkdft
from lattice_models import create_lattice_model


class DMFTPostSolver(DMFTCoreSolver):
    def __init__(self, seedname, params, output_file='', output_group='dmft_out'):

        super(DMFTPostSolver, self).__init__(seedname, params, output_file, output_group, read_only=True, restart=True)


    def calc_dos(self, Sigma_w_sh, mesh, broadening):
        """

        Compute dos in real frequency.

        :param Sigma_w_sh: list
           List of real-frequency self-energy

        :param broadening: float
           Broadening factor

        :return: tuple
           Results are 'dos', 'dosproj', 'dosproj_orb'.

        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'dos'
        params['mu'] = self._chemical_potential
        params['Sigma_w_sh'] = Sigma_w_sh
        params['mesh'] = mesh
        params['broadening'] = broadening
        r = sumkdft.run(os.path.abspath(self._seedname+'.h5'), './work/sumkdft_dos', self._mpirun_command, params)
        return r['dos'], r['dosproj'], r['dosproj_orb']

    def calc_spaghettis(self, Sigma_w_sh, mesh, broadening):
        """

        Compute A(k, omega)

        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'spaghettis'
        params['mu'] = self._chemical_potential
        params['Sigma_w_sh'] = Sigma_w_sh
        params['mesh'] = mesh
        params['broadening'] = broadening
        r = sumkdft.run(os.path.abspath(self._seedname+'.h5'), './work/sumkdft_spaghettis', self._mpirun_command, params)
        return r['akw']

    def calc_momentum_distribution(self):
        """

        Compute momentum distributions and eigen values of H(k)
        Data are taken from bands_data.

        """

        params = self._make_sumkdft_params()
        params['calc_mode'] = 'momentum_distribution'
        params['mu'] = self._chemical_potential
        r = sumkdft.run(os.path.abspath(self._seedname+'.h5'), './work/sumkdft_momentum_distribution', self._mpirun_command, params)
        return r['den']

    def calc_Sigma_w(self, mesh):
        """
        Compute Sigma_w for computing band structure
        For an imaginary-time solver, a list of Nones will be returned.

        :param mesh: (float, float, int)
            real-frequency mesh (min, max, num_points)

        :return: list of Sigma_w

        """

        solver_name = self._params['impurity_solver']['name']
        Solver = impurity_solvers.solver_classes[solver_name]
        if Solver.is_gf_realomega_available():
            Gloc_iw_sh, _ = self.calc_Gloc()
            _, _, sigma_w = self.solve_impurity_models(Gloc_iw_sh, -1, mesh)
            return sigma_w
        else:
            return [None] * self.n_inequiv_shells


class DMFTCoreTools:
    def __init__(self, seedname, params, xk, prefix):
        """
        Class of posting tool for DCore.

        Parameters
        ----------
        :param seedname: string
            name for hdf5 file
        :param params:  dictionary
            Input parameters
        :param xk:  integer array
            x-position for plotting band
        """

        self._params = copy.deepcopy(params)
        # Construct a SumKDFT object
        self._n_pade = int((params['system']['beta']*params['tool']['omega_pade']+numpy.pi) / (2*numpy.pi))
        self._omega_min = float(params['tool']['omega_min'])
        self._omega_max = float(params['tool']['omega_max'])
        self._Nomega = int(params['tool']['Nomega'])
        self._broadening = float(params['tool']['broadening'])
        self._eta = float(params['tool']['eta'])
        self._seedname = seedname
        self._xk = xk
        self._prefix = prefix

        self._solver = DMFTPostSolver(seedname, self._params, output_file=seedname+'.out.h5')
        print("iteration :", self._solver.iteration_number)

    def print_dos(self, dos, dosproj_orb, filename):
        """

        Print DOS to file

        """
        nsh = self._solver.n_inequiv_shells
        om_mesh = numpy.linspace(self._omega_min, self._omega_max, self._Nomega)
        spin_block_names = self._solver.spin_block_names
        inequiv_to_corr = self._solver.inequiv_to_corr
        corr_shell_info = [self._solver.corr_shell_info(ish) for ish in range(self._solver._n_corr_shells)]

        with open(filename, 'w') as f:
            #
            # Description of columns
            #
            print("# [1] Energy", file=f)
            ii = 1
            for isp in spin_block_names:
                ii += 1
                print("# [%d] Total DOS of spin %s" % (ii, isp), file=f)
            for ish in range(nsh):
                block_dim = corr_shell_info[inequiv_to_corr[ish]]['block_dim']
                for isp in spin_block_names:
                    for iorb in range(block_dim):
                        ii += 1
                        print("# [%d] PDOS of shell%d,spin %s,band%d" % (ii, ish, isp, iorb), file=f)
            #
            for iom in range(self._Nomega):
                print("%f" % om_mesh[iom], file=f, end="")
                for isp in spin_block_names:
                    print(" %f" % dos[isp][iom], file=f, end="")
                for ish in range(nsh):
                    block_dim = corr_shell_info[inequiv_to_corr[ish]]['block_dim']
                    for isp in spin_block_names:
                        for iorb in range(block_dim):
                            print(" %f" % dosproj_orb[ish][isp][iom, iorb, iorb].real, end="", file=f)
                print("", file=f)
        print("\n    Output {0}".format(filename))

    def print_band(self, akw, filename):
        """
        Print A(k,w) into a file

        Parameters
        ----------
        akw
        filename

        """

        om_mesh = numpy.linspace(self._omega_min, self._omega_max, self._Nomega)
        with open(filename, 'w') as f:
            offset = 0.0
            for isp in self._solver.spin_block_names:
                # for ik in range(self._n_k):
                for ik, xk in enumerate(self._xk):
                    for iom, om in enumerate(om_mesh):
                        print("%f %f %f" % (xk+offset, om, akw[isp][ik, iom]), file=f)
                    print("", file=f)
                offset = self._xk[-1] * 1.1
                print("", file=f)
        print("\n    Output {0}".format(filename))

    def post(self):
        """
        Calculate DOS (Density Of State) and energy dispersions.
        For Hubbard-I solver, self-energy is calculated in this function.
        For cthyb (both TRIQS and ALPS), self-energy is read from hdf5 file.
        """

        print("\n#############  Compute Green's Function in the Real Frequency  ################\n")

        #
        # Real-frequency self-energy
        #
        mesh = [self._omega_min, self._omega_max, self._Nomega]
        sigma_w_sh = self._solver.calc_Sigma_w(mesh)
        Sigma_iw_sh = self._solver.Sigma_iw_sh(self._solver.iteration_number)
        for ish in range(self._solver.n_inequiv_shells):
            if not sigma_w_sh[ish] is None:
                continue

            # set BlockGf sigma_w
            Sigma_iw = Sigma_iw_sh[ish]
            block_names = self._solver.spin_block_names
            def glist():
                return [GfReFreq(indices=sigma.indices, window=(self._omega_min, self._omega_max),
                                 n_points=self._Nomega, name="sig_pade") for block, sigma in Sigma_iw]
            sigma_w_sh[ish] = BlockGf(name_list=block_names, block_list=glist(), make_copies=False)
            # Analytic continuation
            for bname, sig in Sigma_iw:
                sigma_w_sh[ish][bname].set_from_pade(sig, n_points=self._n_pade, freq_offset=self._eta)

        #
        #  (Partial) DOS
        #
        print("\n#############  Compute (partial) DOS  ################\n")
        dos, dosproj, dosproj_orb = self._solver.calc_dos(sigma_w_sh, mesh, self._broadening)
        self.print_dos(dos, dosproj_orb, self._prefix + self._seedname+'_dos.dat')

        #
        # Band structure
        #
        if self._xk is None:
            return
        #
        print("\n#############  Compute Band Structure  ################\n")
        akw = self._solver.calc_spaghettis(sigma_w_sh, mesh, self._broadening)
        #
        # Print band-structure into file
        #
        self.print_band(akw, self._prefix + self._seedname + '_akw.dat')

    def momentum_distribution(self):
        """
        Calculate Momentum distribution
        """
        if self._xk is None:
            return

        print("\n#############  Momentum Distribution  ################\n")

        den = self._solver.calc_momentum_distribution()

        spn = self._solver.spin_block_names

        n_k, n_orbitals = den.shape[0], den.shape[2]

        SO = 1 if self._solver.use_spin_orbit else 0

        #
        # Output momentum distribution to file
        #
        filename = self._prefix + self._seedname + "_momdist.dat"
        print("\n Output Momentum distribution : ", filename)
        with open(filename, 'w') as fo:
            print("# Momentum distribution", file=fo)
            #
            # Column information
            #
            print("# [Column] Data", file=fo)
            print("# [1] Distance along k-path", file=fo)
            icol = 1
            for isp in spn:
                for iorb in range(n_orbitals):
                    for jorb in range(n_orbitals):
                        icol += 1
                        print("# [%d] Re(MomDist_{spin=%s, %d, %d})" % (icol, isp, iorb, jorb), file=fo)
                        icol += 1
                        print("# [%d] Im(MomDist_{spin=%s, %d, %d})" % (icol, isp, iorb, jorb), file=fo)
            #
            # Write data
            #
            for ik in range(n_k):
                print("%f " % self._xk[ik], end="", file=fo)
                for isp in range(2-SO):
                    for iorb in range(n_orbitals):
                        for jorb in range(n_orbitals):
                            print("%f %f " % (den[ik, isp, iorb, jorb].real,
                                              den[ik, isp, iorb, jorb].imag), end="", file=fo)
                print("", file=fo)


def __print_paramter(p, param_name):
    """
    Print parameters.

    Parameters
    ----------
    p : dictionary
        Dictionary for parameters
    param_name : string
        key for p
    """
    print(param_name + " = " + str(p[param_name]))


def gen_script_gnuplot(xnode, seedname, prefix, spin_orbit):
    file_akw_gp = prefix + seedname + '_akw.gp'

    def print_klabel(label, x, f, with_comma=True):
        print('  "{}"  {}'.format(label, x), end="", file=f)
        if with_comma:
            print(',', end="", file=f)
        print(' \\', file=f)

    k_end = len(xnode) - 1

    with open(file_akw_gp, 'w') as f:
        print("set size 0.95, 1.0", file=f)
        print("set xtics (\\", file=f)
        if spin_orbit:
            for i, node in enumerate(xnode):
                print_klabel(node.label, node.x, f, i != k_end)
        else:
            for node in xnode:
                print_klabel(node.label, node.x, f)
            offset = xnode[-1].x * 1.1
            for i, node in enumerate(xnode):
                print_klabel(node.label, node.x + offset, f, i != k_end)
        print("  )", file=f)
        print("set pm3d map", file=f)
        print("#set pm3d interpolate 5, 5", file=f)
        print("unset key", file=f)
        print("set ylabel \"Energy\"", file=f)
        print("set cblabel \"A(k,w)\"", file=f)
        print("splot \"{0}_akw.dat\"".format(seedname), file=f)
        print("pause -1", file=f)

        print("    Usage:")
        print("\n      $ gnuplot {0}".format(os.path.basename(file_akw_gp)))


def dcore_post(filename, np=1, prefix="./"):
    """
    Main routine for the post-processing tool

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
    pars = create_parser(['model', 'system', 'impurity_solver', 'tool', 'mpi'])
    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    seedname = p["model"]["seedname"]
    p["mpi"]["num_processes"] = np
    mpirun_command = p['mpi']['command'].replace('#', str(p['mpi']['num_processes']))
    mpirun_command_np1 = p['mpi']['command'].replace('#', '1')

    # make directory
    dir = os.path.dirname(prefix)
    if not os.path.exists(dir):
        os.makedirs(dir)

    #
    # Summary of input parameters
    #
    print("\n  @ Parameter summary")
    print("\n    [model] block")
    for k, v in p["model"].items():
        print("      {0} = {1}".format(k, v))
    print("\n    [tool] block")
    for k, v in p["tool"].items():
        print("      {0} = {1}".format(k, v))

    #
    # Generate k-path and compute H(k) on this path
    #
    print("\n################  Generating k-path  ##################\n")

    lattice_model = create_lattice_model(p)
    xk, xnode = lattice_model.generate_Hk_path(p)

    if xk is None:
        print('  A(k,w) calc will be skipped')
    else:
        print("   Total number of k =", len(xk))
        print("    k-point  x")
        for node in xnode:
            print("     %6s  %f" %(node.label, node.x))

    #
    # Coompute DOS and A(k,w)
    #
    print("\n#############   Run DMFTCoreTools  ########################\n")
    dct = DMFTCoreTools(seedname, p, xk, prefix)
    dct.post()
    dct.momentum_distribution()

    #
    # Output gnuplot script
    #
    if xnode is not None:
        print("\n#############   Generate GnuPlot Script  ########################\n")
        gen_script_gnuplot(xnode, seedname, prefix, p["model"]["spin_orbit"])

    print("\n#################  Done  #####################\n")


if __name__ == '__main__':
    from .option_tables import generate_all_description
    import argparse
    from .version import version

    parser = argparse.ArgumentParser(
        prog='dcore_post.py',
        description='pre script for dcore.',
        usage='$ dcore_post input --np 4',
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
    parser.add_argument('--prefix',
                        action='store',
                        default='post/',
                        type=str,
                        help='prefix for output files'
                        )

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)
    dcore_post(args.path_input_file, int(args.np), args.prefix)
