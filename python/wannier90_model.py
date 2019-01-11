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

import argparse
import numpy
import sys
import re

from pytriqs.applications.dft.converters.wannier90_converter import Wannier90Converter
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi

def __generate_wannier90_model(seedname, ncor, norb, n_k, kvec):
    """
    Compute hopping etc. for A(k,w) of Wannier90

    Parameters
    ----------
    seedname : str
    ncor : int
    norb : int
    n_k : integer
        Number of k points
    kvec : float array
        k-points where A(k,w) is computed

    Returns
    -------
    hopping : complex
        k-dependent one-body Hamiltonian
    n_orbitals : integer
        Number of orbitals at each k. It does not depend on k
    proj_mat : complex
        Projection onto each correlated orbitals
    """
    n_spin = 1
    norb_list = re.findall(r'\d+', norb)
    norb = [int(norb_list[icor]) for icor in range(ncor)]
    #
    if mpi.is_master_node():
        print("               ncor = ", ncor)
        for i in range(ncor):
            print("     norb[{0}] = {1}".format(i, norb[i]))
    #
    # Read hopping in the real space from the Wannier90 output
    #
    w90c = Wannier90Converter(seedname=seedname)
    nr, rvec, rdeg, nwan, hamr = w90c.read_wannier90hr(seedname+"_hr.dat")
    #
    # Fourier transformation of the one-body Hamiltonian
    #
    n_orbitals = numpy.ones([n_k, n_spin], numpy.int) * nwan
    hopping = numpy.zeros([n_k, n_spin, numpy.max(n_orbitals), numpy.max(n_orbitals)], numpy.complex_)
    for ik in range(n_k):
        for ir in range(nr):
            rdotk = numpy.dot(kvec[ik, :], rvec[ir, :])
            factor = (numpy.cos(rdotk) + 1j * numpy.sin(rdotk)) / float(rdeg[ir])
            hopping[ik, 0, :, :] += factor * hamr[ir][:, :]
    #
    # proj_mat is (norb*norb) identities at each correlation shell
    #
    proj_mat = numpy.zeros([n_k, n_spin, ncor, numpy.max(norb), numpy.max(n_orbitals)], numpy.complex_)
    iorb = 0
    for icor in range(ncor):
        proj_mat[:, :, icor, 0:norb[icor], iorb:iorb + norb[icor]] = numpy.identity(norb[icor], numpy.complex_)
        iorb += norb[icor]

    return hopping, n_orbitals, proj_mat


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(
            description='Internal program for launching Wannier90Converter')
        parser.add_argument('input_file')
        parser.add_argument('output_file')
        args = parser.parse_args()

        with HDFArchive(args.input_file, 'r') as h:
            seedname = h['seedname']
            ncor = h['ncor']
            norb = h['norb']
            n_k = h['n_k']
            kvec = h['kvec']

        hopping, n_orbitals, proj_mat = __generate_wannier90_model(seedname, ncor, norb, n_k, kvec)

        with HDFArchive(args.output_file, 'w') as h:
            h['hopping'] = hopping
            h['n_orbitals'] = n_orbitals
            h['proj_mat'] = proj_mat

    except Exception as e:
        print("Unexpected error:", e)
        import traceback
        traceback.print_exc()
        sys.exit(1)
