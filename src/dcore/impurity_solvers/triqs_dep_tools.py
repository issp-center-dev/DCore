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

from triqs.operators import *
from itertools import product

def creat_mapping_flatten_index(gf_struct):
    # Map (block_name, index) to an index in the flatten spin-orbital space
    # If blocks are 'up' and 'down', 'up' appears FIRST.
    if isinstance(gf_struct, list):
        gf_struct = {x[0]: x[1] for x in gf_struct}

    to_flatten_index = {}
    if len(gf_struct) == 1:
        from_flatten_index = []
        for offset, index in enumerate(gf_struct['ud']):
            to_flatten_index[('ud', index)] = offset
            from_flatten_index.append(('ud', index))
    else:
        block_names = ['up', 'down']
        from_flatten_index = []
        offset = 0
        for name in block_names:
            for index in gf_struct[name]:
                to_flatten_index[(name, index)] = offset
                from_flatten_index.append((name, index))
                offset += 1

    return to_flatten_index, from_flatten_index


def make_h_int(u_mat, gf_struct):
    """
    Construct an operator representing the interacting Hamiltonian

    :param u_mat: four-index U matrix.
        The dimensions of each axis is spin * orbital.
           gf_struct: dict
    """

    n_orb = int(u_mat.shape[0]/2)
    _, from_flatten_index = creat_mapping_flatten_index(gf_struct)

    ham = Operator()
    for i1, i2, i3, i4 in product(list(range(2*n_orb)), repeat=4):
        ham += 0.5 * u_mat[i1, i2, i3, i4] \
               * c_dag(*from_flatten_index[i1]) * c_dag(*from_flatten_index[i2]) \
               * c(*from_flatten_index[i4]) * c(*from_flatten_index[i3])

    return ham