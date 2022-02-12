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
import sys, os
from dcore._dispatcher import MeshReFreq, MeshImFreq, GfReFreq, GfImFreq
from dcore.program_options import create_parser, parse_parameters
from dcore.version import version, print_header
from spm_omega import AnaContSmooth, AnaContSpM
import numpy
import toml


def dcore_spm(filename):
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
    pars = create_parser(['tool'])

    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    seedname = p["model"]["seedname"]
    niter = p["tool"]["niter_spm"]
    alpha_min = p["tool"]["alpha_min_spm"]
    alpha_max = p["tool"]["alpha_max_spm"]


    print("Reading ", seedname + "_anacont.toml...")
    with open(seedname + "_anacont.toml", "r") as f:
        params = toml.load(f)
    beta = params['beta']

    print("Reading ", seedname + "_sigma_iw.npz...")
    npz = numpy.load(seedname + "_sigma_iw.npz")

    # Real-frequency mesh
    assert params["omega_min"] < params["omega_max"]
    wmax = max(numpy.abs(params["omega_min"]), numpy.abs(params["omega_max"]))
    mesh_w = numpy.linspace(params["omega_min"], params["omega_max"], params["Nomega"])

    data_w = {}

    num_data = numpy.sum([key.startswith("data") for key in npz.keys()])
    for idata in range(num_data):
        key = f"data{idata}"
        key_hf = f"hartree_fock{idata}"
        sigma_iw = npz[key]
        hf = npz[key_hf]
        nw = sigma_iw.shape[0]//2

        # Remove Hartree-Fock term
        sigma_iw -= hf[:, None, None]

        vsample = 2*numpy.arange(-nw, nw)+1
        solver = AnaContSmooth(beta, wmax, "F", "freq", vsample)
        x, _ = solver.solve_elbow(sigma_iw, alpha_min, alpha_max, niter=niter)

        data_w[key] = solver.g_omega(x, mesh_w) + hf[:, None, None]

    print("Writing to", seedname + "_sigma_w.npz...")
    numpy.savez(seedname + "_sigma_w.npz", **data_w)


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version, print_header

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_spm.py',
        description='Analytic continuation of self-energy using sparse-modeling techniques',
        usage='$ dcore_spm input',
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
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)
    dcore_spm(args.path_input_file)