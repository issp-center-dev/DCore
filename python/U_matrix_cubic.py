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
from pytriqs.operators.util.U_matrix import *

def U_matrix_kanamori_4index_t2g(U_int, J_hund):
    F = [U_int, 63.*J_hund, -63.*J_hund]
    Umat = U_matrix(l=2, radial_integrals=F, basis='cubic')
    U_sub = t2g_submatrix(Umat)
    return U_sub

def U_matrix_kanamori_4index_eg(U_int, J_hund):
    F = [U_int, 21.*J_hund, -21.*J_hund]
    Umat = U_matrix(l=2, radial_integrals=F, basis='cubic')
    U_sub = eg_submatrix(Umat)
    return U_sub

def U_matrix_kanamori_4index_p(U_int, J_hund):
    F = [U_int-4./3.*J_hund, 25./3.*J_hund]
    Umat = U_matrix(l=1, radial_integrals=F, basis='cubic')
    return Umat
