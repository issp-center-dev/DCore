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


import numpy
import h5py

from dcore.dc import *

def test_dc():
    # U n_up n_down
    U = 1.0
    nso = 2
    u_mat = numpy.zeros((nso, nso, nso, nso), dtype=numpy.complex128)
    u_mat[0, 1, 0, 1] = 2*U

    dm = numpy.zeros((nso, nso), dtype=numpy.complex128)
    # Only up spin
    dm[0,0] = 1.0

    numpy.testing.assert_allclose(hartree_fock_term(dm, u_mat), numpy.diag([0, 1]))


def test_dc_dict():
    n_up = 0.1
    n_dn = 0.2
    U = 1.0
    nso = 2
    u_mat = numpy.zeros((nso, nso, nso, nso), dtype=numpy.complex128)
    u_mat[0, 1, 0, 1] = 2*U

    dm = {"up": n_up * numpy.identity(1), "down": n_dn * numpy.identity(1)}
    dc = hf_dc(dm, u_mat, False)
    numpy.testing.assert_allclose(dc["up"], U*n_dn*numpy.identity(1))
    numpy.testing.assert_allclose(dc["down"], U*n_up*numpy.identity(1))

    dm_so = {"ud": numpy.diag([n_up, n_dn])}
    dc_so = hf_dc(dm_so, u_mat, True)
    numpy.testing.assert_allclose(dc_so["ud"], numpy.diag([U*n_dn, U*n_up]))
