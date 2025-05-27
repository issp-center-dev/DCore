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

import argparse
import sys
import os.path
import itertools
import numpy
from dcore._dispatcher import MeshReFreq
from dcore.version import version, print_header
from dcore.program_options import create_parser, parse_parameters
from dcore.anacont import pade, spm


def dcore_anacont(inifile):
    # tool is removed but read for error message
    parser = create_parser(["system", "model", "tool", "post", "post.anacont"])
    parser.read(inifile)
    params = parser.as_dict()
    parse_parameters(params)

    beta = params["system"]["beta"]

    omega_min = params["post.anacont"]["omega_min"]
    omega_max = params["post.anacont"]["omega_max"]
    Nomega = params["post.anacont"]["Nomega"]
    mesh_w = MeshReFreq(omega_min, omega_max, Nomega)

    seedname = params["model"]["seedname"]
    file_sigma_iw = seedname + "_sigma_iw.npz"
    # file_sigma_iw = params["post"]["file_sigma_iw"]
    if not os.path.exists(file_sigma_iw):
        sys.exit("File not found: " + file_sigma_iw)
    sigma_iw_npz = numpy.load(file_sigma_iw)

    dir_work = params["post"]["dir_work"]
    if not os.path.exists(dir_work):
        os.makedirs(dir_work)

    solver = params["post.anacont"]["solver"]
    if solver == "pade":
        params_ac = pade.parameters_from_ini(inifile)
        data_w = pade.anacont(sigma_iw_npz, beta, mesh_w, params_ac)
    elif solver == "spm":
        params_ac = spm.parameters_from_ini(inifile)
        if "post.anacont.spm.solver" not in params_ac:
            params_ac["post.anacont.spm.solver"] = {}
        data_w = spm.anacont(
            sigma_iw_npz,
            beta,
            mesh_w,
            params_ac["post.anacont.spm"],
            params_ac["post.anacont.spm.solver"],
        )
    else:
        sys.exit("ERROR: Unknown anacont solver " + solver)
    omega = numpy.array(list(mesh_w.values()))
    data_w["omega"] = omega

    dir_post = params["post"]["dir_post"]
    if not os.path.exists(dir_post):
        os.makedirs(dir_post)
    file_sigma_w = os.path.join(dir_post, "sigma_w.npz")
    print("Writing to", file_sigma_w + "...")
    numpy.savez(file_sigma_w, **data_w)

    if params["post.anacont"]["show_result"] or params["post.anacont"]["save_result"]:
        import matplotlib.pyplot as plt

        ndata = 0
        for key in data_w.keys():
            if key.startswith("data"):
                ndata += 1
        for idata in range(ndata):
            sigma_w = data_w[f"data{idata}"]
            for iorb, jorb in itertools.product(
                range(sigma_w.shape[1]), range(sigma_w.shape[2])
            ):
                plt.axhline(y=0, xmin=omega_min, xmax=omega_max, color="lightgrey")
                plt.plot(omega, sigma_w[:, iorb, jorb].real, label="Real")
                plt.plot(omega, sigma_w[:, iorb, jorb].imag, label="Imag")
                plt.xlim(omega_min, omega_max)
                plt.xlabel(r"$\omega$")
                plt.legend()
                plt.title(rf"$\Sigma_{{{iorb}{jorb}}}( \omega )$ of shell {idata}")
                plt.tight_layout()
                if params["post.anacont"]["save_result"]:
                    file_result = os.path.join(
                        dir_work, f"sigma_w_{idata}_{iorb}_{jorb}.png"
                    )
                    print("Writing to", file_result + "...")
                    plt.savefig(file_result)
                if params["post.anacont"]["show_result"]:
                    plt.show()
                plt.close()


def run():

    print_header()

    parser = argparse.ArgumentParser(
        prog="dcore_anacont.py",
        description="DCore post script -- analytic continuation.",
        usage="$ dcore_anacont input.ini",
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        # epilog=generate_all_description()
    )
    parser.add_argument(
        "path_input_files",
        action="store",
        default=None,
        type=str,
        nargs="*",
        help="Input filename(s)",
    )
    parser.add_argument(
        "--version", action="version", version="DCore {}".format(version)
    )

    args = parser.parse_args()

    dcore_anacont(args.path_input_files)
