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

import copy
import argparse
import itertools
from dcore._dispatcher import MeshReFreq, MeshImFreq, GfReFreq, GfImFreq
from dcore.fourier import _fft_fermion_w2t, _matsubara_freq_fermion
from dcore.version import version, print_header
from dcore.tools import calc_retarted_function_from_spectrum
import numpy
import toml

import spmac.main as spmac


def dcore_pyspmac(seedname):
    print("Reading ", seedname + "_anacont.toml...")
    with open(seedname + "_anacont.toml", "r") as f:
        params = toml.load(f)

    print("Reading ", seedname + "_sigma_iw.npz...")
    npz = numpy.load(seedname + "_sigma_iw.npz")

    assert params["omega_min"] < params["omega_max"]
    # mesh_w = MeshReFreq(params["omega_min"], params["omega_max"], params["Nomega"])

    beta = params["beta"]
    omega_min = params["omega_min"]
    omega_max = params["omega_max"]
    Nomega = params["Nomega"]
    ws_rho = numpy.linspace(omega_min, omega_max, num=Nomega)
    ws_sw = copy.copy(ws_rho)
    ws_rho[0] -= 1e-8
    ws_rho[-1] += 1e-8

    params_spmac = params.get("pyspmac", {})
    params_spmac["min_omega"] = omega_min
    params_spmac["max_omega"] = omega_max
    params_spmac["num_omega"] = Nomega
    params_spmac["beta"] = beta

    data_w = {}

    num_data = numpy.sum([key.startswith("data") for key in npz.keys()])
    for idata in range(num_data):
        key = f"data{idata}"
        siw = npz[key]
        hf: complex = npz[f"hartree_fock{idata}"]
        siw -= hf
        nt, nflavors, _ = siw.shape
        niw = nt // 2
        inv_iw = 1.0 / numpy.imag(_matsubara_freq_fermion(beta, niw))
        ntail = 4

        params_spmac["num_flavor"] = nflavors
        st = numpy.zeros((0, 0, 0), dtype=numpy.float64)
        for i, j in itertools.product(range(nflavors), repeat=2):
            s = siw[:, i, j]
            ax = -numpy.sum(numpy.imag(s[-ntail:]) * inv_iw[-ntail:]) / numpy.sum(
                inv_iw[-ntail:] ** 2
            )
            st_ij = _fft_fermion_w2t(siw[:, i, j], beta, a=ax)
            if st.shape[2] == 0:
                st = numpy.zeros((nflavors, nflavors, st_ij.size), dtype=numpy.float64)
            st[i, j, :] = numpy.real(st_ij[:])
        params_spmac["output"] = f"output{idata}"
        solver, rho, loglambda = spmac.run(params_spmac, Gtau=st)
        solver.write_rho(solver.ws, rho, loglambda=loglambda, outdir=solver.outdir)
        if nflavors == 1:
            rho = rho.reshape(-1, 1, 1)
        else:
            rho = rho.transpose(1, 2, 0)
        G_w = calc_retarted_function_from_spectrum(
            ws_sw, ws_rho, rho, rho_is_complex=False
        )
        data_w[key] = G_w + hf
    print("Writing to", seedname + "_sigma_w.npz...")
    numpy.savez(seedname + "_sigma_w.npz", **data_w)


def run():

    print_header()

    parser = argparse.ArgumentParser(
        prog="dcore_pyspmac.py",
        description="pre script for dcore.",
        usage="$ dcore_pyspmac input",
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

    dcore_pyspmac(args.seedname)
