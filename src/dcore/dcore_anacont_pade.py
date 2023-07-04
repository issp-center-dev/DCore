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
from dcore._dispatcher import MeshReFreq, MeshImFreq, GfReFreq, GfImFreq
from dcore.version import version, print_header
import numpy
import toml

def _set_n_pade(omega_cutoff, beta, n_min, n_max):
    """
    Return (int)n_pade: the number of Matsubara frequencies below the cutoff frequency.
    n_pade is bounded between n_min and n_max
    """
    n_pade = int((beta * omega_cutoff + numpy.pi) / (2.0 * numpy.pi))
    print("n_pade = {} (evaluated from omega_pade)".format(n_pade))
    n_pade = max(n_pade, n_min)
    n_pade = min(n_pade, n_max)
    print("n_pade = {}".format(n_pade))
    return n_pade

def dcore_pade(seedname):
    print("Reading ", seedname + "_anacont.toml...")
    with open(seedname + "_anacont.toml", "r") as f:
        params = toml.load(f)
 
    print("Reading ", seedname + "_sigma_iw.npz...")
    npz = numpy.load(seedname + "_sigma_iw.npz")

    assert params["omega_min"] < params["omega_max"]
    mesh_w = MeshReFreq(params["omega_min"], params["omega_max"], params["Nomega"])

    n_pade = _set_n_pade(params["pade"]["omega_max"], params['beta'], n_min=params["pade"]["n_min"], n_max=params["pade"]["n_max"])

    data_w = {}

    num_data = numpy.sum([key.startswith("data") for key in npz.keys()])
    for idata in range(num_data):
        key = f"data{idata}"
        data = npz[key]
        mesh_iw = MeshImFreq(params["beta"], "Fermion", data.shape[0]//2)
        sigma_iw = GfImFreq(data=data, beta=params["beta"], mesh=mesh_iw)
        sigma_w = GfReFreq(mesh=mesh_w, target_shape=data.shape[1:])
        sigma_w.set_from_pade(sigma_iw, n_points=n_pade, freq_offset=params["pade"]["eta"])
        data_w[key] = sigma_w.data
    print("Writing to", seedname + "_sigma_w.npz...")
    numpy.savez(seedname + "_sigma_w.npz", **data_w)

def run():

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_pade.py',
        description='pre script for dcore.',
        usage='$ dcore_pade input',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        #epilog=generate_all_description()
    )
    parser.add_argument('seedname',
                        action='store',
                        default=None,
                        type=str,
                        help="seedname"
                        )
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))

    args = parser.parse_args()

    dcore_pade(args.seedname)
