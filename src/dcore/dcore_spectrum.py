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
import copy
from itertools import product

from dcore._dispatcher import *
from dcore.dmft_core import DMFTCoreSolver
from dcore.program_options import (
    create_parser,
    parse_parameters,
    parse_bvec,
    _set_nk,
    print_parameters,
    delete_parameters,
)
from dcore.tools import save_Sigma_w_sh_txt
from dcore import impurity_solvers
from dcore.lattice_models import create_lattice_model
from .sumkdft_workers.launcher import run_sumkdft

def _read_ws(npz_file):
    npz = np.load(npz_file)
    ws = npz["omega"]
    return ws

def _read_sigma_w(npz_file, nsh, mesh, block_names):
    npz = np.load(npz_file)
    Sigma_w = []
    idx = 0
    for _ in range(nsh):
        block_list = []
        for _ in range(len(block_names)):
            block_list.append(GfReFreq(data=npz[f"data{idx}"], mesh=mesh))
            idx += 1
        G = BlockGf(name_list=block_names, block_list=block_list, make_copies=False)
        Sigma_w.append(G)
    return Sigma_w


class DMFTPostSolver(DMFTCoreSolver):
    def __init__(self, seedname, params, output_file="", output_group="dmft_out"):

        super(DMFTPostSolver, self).__init__(
            seedname, params, output_file, output_group, read_only=True, restart=True
        )

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
        params["mu"] = self._chemical_potential
        params["Sigma_w_sh"] = Sigma_w_sh
        params["mesh"] = mesh
        params["broadening"] = broadening
        r = run_sumkdft(
            "SumkDFTWorkerDOS",
            os.path.abspath(self._seedname + ".h5"),
            "./work/sumkdft_dos",
            self._mpirun_command,
            params,
        )
        return r["dos"], r["dosproj"], r["dosproj_orb"]

    def calc_spaghettis(self, Sigma_w_sh, mesh, broadening, kmesh_type):
        """

        Compute A(k, omega)

        """

        params = self._make_sumkdft_params()
        params["calc_mode"] = "spaghettis"
        params["mu"] = self._chemical_potential
        params["Sigma_w_sh"] = Sigma_w_sh
        params["mesh"] = mesh
        params["broadening"] = broadening
        if kmesh_type == "line":
            params["bands_data"] = "dft_bands_input"
        elif kmesh_type == "mesh":
            params["bands_data"] = "dft_bands_mesh_input"
        else:
            raise RuntimeError("Invalid kmesh_type: {}".format(kmesh_type))
        r = run_sumkdft(
            "SumkDFTWorkerSpaghettis",
            os.path.abspath(self._seedname + ".h5"),
            "./work/sumkdft_spaghettis",
            self._mpirun_command,
            params,
        )
        return r["akw"]

    def calc_momentum_distribution(self):
        """

        Compute momentum distributions and eigen values of H(k)
        Data are taken from bands_data.

        """

        params = self._make_sumkdft_params()
        params["calc_mode"] = "momentum_distribution"
        params["mu"] = self._chemical_potential
        r = run_sumkdft(
            "SumkDFTWorkerMomDist",
            os.path.abspath(self._seedname + ".h5"),
            "./work/sumkdft_mom_dist",
            self._mpirun_command,
            params,
        )
        return r["den"]


class DMFTCoreTools:
    def __init__(self, seedname, params, n_k, xk, nkdiv_mesh, kvec_mesh, dir_post):
        """
        Class of posting tool for DCore.

        Parameters
        ----------
        :param seedname: string
            name for hdf5 file
        :param params:  dictionary
            Input parameters
        :param n_k: integer
            Number of k points
        :param xk:  integer array
            x-position for plotting band
        :param nkdiv_mesh:  (int, int, int)
            Number of k points along each axis for computing A(k, omega)
        :param kvec_mesh:  float array of dimension (*, 3)
            k points in fractional coordinates for computing A(k, omega) on a mesh
        """

        self._params = copy.deepcopy(params)
        self._sigma_w_file = os.path.join(dir_post, "sigma_w.npz")
        ws = _read_ws(self._sigma_w_file)
        self.mesh_w = MeshReFreq(ws[0], ws[-1], len(ws))
        if TRIQS_COMPAT:
            self.mesh_w._points = ws
        else:
            ws_linear = np.linspace(ws[0], ws[-1], num=len(ws))
            if not np.allclose(ws, ws_linear):
                sys.exit("Error: omega points must be linearly spaced when using TRIQS (DCORE_TRIQS_COMPAT=0)")

        # Construct a SumKDFT object
        self._omega_min = ws[0]
        self._omega_max = ws[-1]
        self._Nomega = len(ws)
        # self._omega_min = float(params["post"]["omega_min"])
        # self._omega_max = float(params["post"]["omega_max"])
        # self._Nomega = int(params["post"]["Nomega"])
        self._broadening = float(params["post.spectrum"]["broadening"])
        self._seedname = seedname
        self._n_k = n_k
        self._xk = xk
        self._kvec_mesh = kvec_mesh
        self._nkdiv_mesh = nkdiv_mesh
        self._dir_post = dir_post


        self._solver = DMFTPostSolver(
            seedname, self._params, output_file=seedname + ".out.h5"
        )
        print("iteration :", self._solver.iteration_number)

    def print_dos(self, dos, dosproj_orb, filename):
        """

        Print DOS to file

        """
        nsh = self._solver.n_inequiv_shells
        om_mesh = numpy.linspace(self._omega_min, self._omega_max, self._Nomega)
        spin_block_names = self._solver.spin_block_names
        inequiv_to_corr = self._solver.inequiv_to_corr
        corr_shell_info = [
            self._solver.corr_shell_info(ish)
            for ish in range(self._solver._n_corr_shells)
        ]

        filepath = os.path.join(self._dir_post, filename)
        with open(filepath, "w") as f:
            #
            # Description of columns
            #
            print("# [1] Energy", file=f)
            ii = 1
            for isp in spin_block_names:
                ii += 1
                print("# [%d] Total DOS of spin %s" % (ii, isp), file=f)
            for ish in range(nsh):
                block_dim = corr_shell_info[inequiv_to_corr[ish]]["block_dim"]
                for isp in spin_block_names:
                    for iorb in range(block_dim):
                        ii += 1
                        print(
                            "# [%d] PDOS of shell%d,spin %s,band%d"
                            % (ii, ish, isp, iorb),
                            file=f,
                        )
            #
            for iom in range(self._Nomega):
                print("%f" % om_mesh[iom], file=f, end="")
                for isp in spin_block_names:
                    print(" %f" % dos[isp][iom], file=f, end="")
                for ish in range(nsh):
                    block_dim = corr_shell_info[inequiv_to_corr[ish]]["block_dim"]
                    for isp in spin_block_names:
                        for iorb in range(block_dim):
                            print(
                                " %f" % dosproj_orb[ish][isp][iom, iorb, iorb].real,
                                end="",
                                file=f,
                            )
                print("", file=f)
        print("\n    Output {0}".format(filepath))

    def print_band(self, akw, filename):
        """
        Print A(k,w) into a file

        Parameters
        ----------
        akw
        filename

        """

        filepath = os.path.join(self._dir_post, filename)
        om_mesh = numpy.linspace(self._omega_min, self._omega_max, self._Nomega)
        with open(filepath, "w") as f:
            offset = 0.0
            for isp in self._solver.spin_block_names:
                for ik, xk in enumerate(self._xk):
                    for iom, om in enumerate(om_mesh):
                        print("%f %f %f" % (xk + offset, om, akw[isp][ik, iom]), file=f)
                    print("", file=f)
                offset = self._xk[-1] * 1.1
                print("", file=f)
        print("\n    Output {0}".format(filepath))

    def post(self):
        """
        Calculate DOS (Density Of State) and energy dispersions.
        For Hubbard-I solver, self-energy is calculated in this function.
        For cthyb (both TRIQS and ALPS), self-energy is read from hdf5 file.
        """

        print(
            "\n#############  Load Green's Function in the Real Frequency  ################\n"
        )

        #
        # Real-frequency self-energy
        #
        if not os.path.exists(self._sigma_w_file):
            sys.exit("File not found: " + self._sigma_w_file)
        sigma_w_sh = _read_sigma_w(
            self._sigma_w_file,
            self._solver.n_inequiv_shells,
            self.mesh_w,
            self._solver.spin_block_names,
        )

        print(
            "\n#############  Print Self energy in the Real Frequency  ################\n"
        )
        filename = os.path.join(self._dir_post, "sigma_w.dat")
        print("\n Writing real-freqnecy self-energy into ", filename)
        save_Sigma_w_sh_txt(filename, sigma_w_sh, self._solver.spin_block_names)

        mesh = [self._omega_min, self._omega_max, self._Nomega]
        #
        #  (Partial) DOS
        #
        print("\n#############  Compute (partial) DOS  ################\n")
        dos, dosproj, dosproj_orb = self._solver.calc_dos(
            sigma_w_sh, mesh, self._broadening
        )
        self.print_dos(dos, dosproj_orb, "dos.dat")

        #
        # Band structure
        #
        if self._xk is not None:
            print("\n#############  Compute Band Structure  ################\n")
            akw = self._solver.calc_spaghettis(
                sigma_w_sh, mesh, self._broadening, "line"
            )
            self.print_band(akw, "akw.dat")

        #
        # A(k, omega) on a mesh
        #
        nk_mesh = numpy.prod(self._nkdiv_mesh)
        if nk_mesh > 0:
            print("\n#############  Compute A(k, omega) on a mesh ################\n")
            akw = self._solver.calc_spaghettis(
                sigma_w_sh, mesh, self._broadening, "mesh"
            )
            # print("debug", mesh)
            # print("debugB", len(mesh))
            # print("debugC", self._Nomega)
            om_mesh = numpy.linspace(mesh[0], mesh[1], mesh[2])
            bvec = parse_bvec(self._params["model"]["bvec"])
            for bname in self._solver.spin_block_names:
                filename = os.path.join(self._dir_post, f"akw_mesh_{bname}.dat")
                print(f"\n    Output {filename}")
                with open(filename, "w") as f:
                    print(
                        "# {} {} {}   {}".format(
                            self._nkdiv_mesh[0],
                            self._nkdiv_mesh[1],
                            self._nkdiv_mesh[2],
                            self._Nomega,
                        ),
                        file=f,
                    )
                    for i in range(3):
                        print("# {} {} {}".format(*bvec[i, :]), file=f)
                    for ik in range(nk_mesh):
                        for iom in range(self._Nomega):
                            print(
                                "%f %f %f    %f %f"
                                % (
                                    self._kvec_mesh[ik, 0],
                                    self._kvec_mesh[ik, 1],
                                    self._kvec_mesh[ik, 2],
                                    om_mesh[iom],
                                    akw[bname][ik, iom],
                                ),
                                file=f,
                            )

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
        filename = os.path.join(self._dir_post, "momdist.dat")
        print("\n Output Momentum distribution : ", filename)
        with open(filename, "w") as fo:
            print("# Momentum distribution", file=fo)
            #
            # Column information
            #
            print("# [Column] Data", file=fo)
            print("# [1] Distance along k-path", file=fo)
            icol = 1

            for isp, iorb, jorb in product(spn, range(n_orbitals), range(n_orbitals)):
                icol += 1
                print(
                    "# [%d] Re(MomDist_{spin=%s, %d, %d})"
                    % (icol, isp, iorb, jorb),
                    file=fo,
                )
                icol += 1
                print(
                    "# [%d] Im(MomDist_{spin=%s, %d, %d})"
                    % (icol, isp, iorb, jorb),
                    file=fo,
                )

            #
            # Write data
            #
            for ik in range(n_k):
                print("%f " % self._xk[ik], end="", file=fo)
                for isp, iorb, jorb in product(range(2-SO), range(n_orbitals), range(n_orbitals)):
                    print(
                        "%f %f "
                        % (
                            den[ik, isp, iorb, jorb].real,
                            den[ik, isp, iorb, jorb].imag,
                        ),
                        end="",
                        file=fo,
                    )
                print("", file=fo)


def __print_parameter(p, param_name):
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


def gen_script_gnuplot(xnode, seedname, dir_post, spin_orbit):
    file_akw_gp = os.path.join(dir_post, "akw.gp")

    def print_klabel(label, x, f, with_comma=True):
        print('  "{}"  {}'.format(label, x), end="", file=f)
        if with_comma:
            print(",", end="", file=f)
        print(" \\", file=f)

    k_end = len(xnode) - 1

    with open(file_akw_gp, "w") as f:
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
        print('set ylabel "Energy"', file=f)
        print('set cblabel "A(k,w)"', file=f)
        print('splot "akw.dat" u 1:2:(abs($3))', file=f)
        print("pause -1", file=f)

        print("    Usage:")
        print("\n      $ gnuplot {0}".format(os.path.basename(file_akw_gp)))


def dcore_spectrum(filename, np=1):
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
    # tool is removed but read for error message
    #
    pars = create_parser(["model", "system", "impurity_solver", "tool", "post", "post.spectrum", "mpi"])
    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)

    seedname = p["model"]["seedname"]
    p["mpi"]["num_processes"] = np
    mpirun_command = p["mpi"]["command"].replace("#", str(p["mpi"]["num_processes"]))
    mpirun_command_np1 = p["mpi"]["command"].replace("#", "1")

    #
    # Delete unnecessary parameters
    #
    delete_parameters(
        p,
        block="model",
        delete=[
            "interaction",
            "density_density",
            "kanamori",
            "slater_f",
            "slater_uj",
            "slater_basis",
            "local_potential_matrix",
            "local_potential_factor",
        ],
    )

    # Summary of input parameters
    print_parameters(p)

    # make directory
    dir_post = p["post"]["dir_post"]


    if dir_post:
        os.makedirs(dir_post, exist_ok=True)
    else:
        dir_post = "."


    #
    # Generate k-path and compute H(k) on this path
    #
    print("\n################  Generating k-path  ##################\n")

    lattice_model = create_lattice_model(p)
    xk, xnode = lattice_model.generate_Hk_path(p)

    if xk is None:
        n_k = 0
        print("  A(k,w) calculation will be skipped")
    else:
        n_k = len(xk)
        print("   Total number of k =", n_k)
        print("    k-point  x")
        for node in xnode:
            print("     %6s  %f" % (node.label, node.x))

    #
    # Generate k mesh and compute H(k) on the mesh
    #
    nk_div = _set_nk(
        p["post.spectrum"]["nk_mesh"],
        p["post.spectrum"]["nk0_mesh"],
        p["post.spectrum"]["nk1_mesh"],
        p["post.spectrum"]["nk2_mesh"],
    )
    kvec_mesh = None
    if all(div != 0 for div in nk_div):
        print(
            "\n################  Constructing H(k) for compute A(k, omega) on a mesh  ##################"
        )
        k = [numpy.linspace(0, 2 * numpy.pi, div + 1)[:-1] for div in nk_div]
        kvec_mesh = numpy.array([kxyz for kxyz in product(k[0], k[1], k[2])])
        lattice_model.write_dft_band_input_data(
            p, kvec_mesh, bands_data="dft_bands_mesh_input"
        )

    #
    # Compute DOS and A(k,w)
    #
    print("\n#############   Run DMFTCoreTools  ########################\n")
    dct = DMFTCoreTools(seedname, p, n_k, xk, nk_div, kvec_mesh, dir_post)
    dct.post()
    dct.momentum_distribution()

    #
    # Output gnuplot script
    #
    if xnode is not None:
        print("\n#############   Generate Gnuplot Script  ########################\n")
        gen_script_gnuplot(xnode, seedname, dir_post, p["model"]["spin_orbit"])

    print("\n#################  Done  #####################\n")


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version, print_header

    print_header()

    parser = argparse.ArgumentParser(
        prog="dcore_spectrum.py",
        description="Post-processing script in DCore",
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description(),
    )
    parser.add_argument(
        "path_input_files",
        action="store",
        default=None,
        type=str,
        nargs="*",
        help="Input filename(s)",
    )
    parser.add_argument("--np", help="Number of MPI processes", required=True)
    parser.add_argument(
        "--version", action="version", version="DCore {}".format(version)
    )

    args = parser.parse_args()

    for path_input_file in args.path_input_files:
        if os.path.isfile(path_input_file) is False:
            sys.exit(f"Input file '{path_input_file}' does not exist.")
    dcore_spectrum(args.path_input_files, int(args.np))
