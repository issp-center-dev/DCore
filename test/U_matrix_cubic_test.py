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

from pytriqs.operators.util.U_matrix import *
from pytriqs.operators.util.hamiltonians import *
# from pytriqs.operators.util.extractors import *
from pytriqs.applications.dcore.U_matrix_cubic import *


def kanamori_4index():
    l=2
    U_int=4.8673 # random
    J_hund=0.8598

    spin_names = ["up", "down"]

    # eg
    orb_names = cubic_names('eg')

    Umat = U_matrix_kanamori_4index_eg(U_int=U_int, J_hund=J_hund)
    h_kanamori_from_slater = h_int_slater(spin_names, orb_names, U_matrix=Umat, off_diag=True)

    U2, U2prime = U_matrix_kanamori(n_orb=2, U_int=U_int, J_hund=J_hund)
    h_kanamori = h_int_kanamori(spin_names, orb_names, U2, U2prime, J_hund, off_diag=True)

    assert (h_kanamori_from_slater - h_kanamori).is_zero()

    # t2g
    orb_names = cubic_names('t2g')

    Umat = U_matrix_kanamori_4index_t2g(U_int=U_int, J_hund=J_hund)
    h_kanamori_from_slater = h_int_slater(spin_names, orb_names, U_matrix=Umat, off_diag=True)

    U2, U2prime = U_matrix_kanamori(n_orb=3, U_int=U_int, J_hund=J_hund)
    h_kanamori = h_int_kanamori(spin_names, orb_names, U2, U2prime, J_hund, off_diag=True)

    assert (h_kanamori_from_slater - h_kanamori).is_zero()

kanamori_4index()
