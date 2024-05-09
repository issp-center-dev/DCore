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
import toml
import numpy as np
from dcore._dispatcher import MeshReFreq
from dcore.version import version, print_header
from dcore.anacont.spm import set_default_values, anacont


def set_default_config(params):
    default_values = {
        "show_fit": False,
        "show_result": False,
        "show_fit": False,
        "verbose_opt": False,
        "max_iters_opt": 100,
        "solver_opt": "ECOS",
    }
    params["spm"] = set_default_values(params["spm"], default_values)
    return params


def dcore_anacont_spm_from_seedname(seedname):
    print("Reading ", seedname + "_anacont.toml...")
    with open(seedname + "_anacont.toml", "r") as f:
        params = toml.load(f)
    params = set_default_config(params)
    print("Using configuration: ", params)

    print("Reading ", seedname + "_sigma_iw.npz...")
    npz = np.load(seedname + "_sigma_iw.npz")

    assert params["omega_min"] < params["omega_max"]
    mesh_w = MeshReFreq(params["omega_min"], params["omega_max"], params["Nomega"])

    assert params["spm"]["n_matsubara"] > 0

    beta = params["beta"]
    params_spm = params["spm"]
    data_w = anacont(npz, beta, mesh_w, params_spm)

    print("Writing to", seedname + "_sigma_w.npz...")
    np.savez(seedname + "_sigma_w.npz", **data_w)


# example file for 'seedname_anacont.toml'
# beta = 40.0
# Nomega = 4000
# omega_min = -6.0
# omega_max = 6.2

# [spm]
# n_matsubara = 300 #number of retained Matsubara frequencies
# n_tail = 100
# n_tau = 10000
# n_sv = 30 # number of retained singular values
# show_fit = false
# show_result = false
# lambda = 1e-6
# verbose_opt = true
# max_iters_opt = 100
# solver_opt = 'ECOS'


def run():
    print_header()

    parser = argparse.ArgumentParser(
        prog="dcore_anacont_spm.py",
        description="pre script for dcore.",
        usage="$ dcore_anacont_spm input",
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

    dcore_anacont_spm_from_seedname(args.seedname)
