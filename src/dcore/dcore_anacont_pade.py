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
import numpy
import toml

from dcore._dispatcher import MeshReFreq
from dcore.version import version, print_header
from dcore.anacont.pade import anacont


def dcore_anacont_pade_from_seedname(seedname):
    print("Reading ", seedname + "_anacont.toml...")
    with open(seedname + "_anacont.toml", "r") as f:
        params = toml.load(f)

    print("Reading ", seedname + "_sigma_iw.npz...")
    sigma_iw_npz = numpy.load(seedname + "_sigma_iw.npz")

    assert params["omega_min"] < params["omega_max"]
    mesh_w = MeshReFreq(params["omega_min"], params["omega_max"], params["Nomega"])

    beta = params["beta"]
    params_pade = params.get("pade", {})
    params_pade["iomega_max"] = params_pade.get("omega_max", -1.0)
    params_pade["eta"] = params_pade.get("eta", 0.01)
    params_pade["n_min"] = params_pade.get("n_min", 0)
    params_pade["n_max"] = params_pade.get("n_max", 100000000)

    data_w = anacont(
        sigma_iw_npz, beta=beta, mesh_w=mesh_w, params_pade=params_pade
    )

    print("Writing to", seedname + "_sigma_w.npz...")
    numpy.savez(seedname + "_sigma_w.npz", **data_w)


def run():

    print_header()

    parser = argparse.ArgumentParser(
        prog="dcore_anacont_pade.py",
        description="pre script for dcore.",
        usage="$ dcore_anacont_pade input",
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        # epilog=generate_all_description()
    )
    parser.add_argument(
        "seedname", action="store", default=None, type=str, help="seedname"
    )
    parser.add_argument(
        "--version", action="version", version="DCore {}".format(version)
    )

    args = parser.parse_args()

    dcore_anacont_pade_from_seedname(args.seedname)
