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
from itertools import *

from typing import Dict

def hartree_fock_term(dm: numpy.ndarray, u_mat: numpy.ndarray) -> numpy.ndarray:
    """Compute Hartree-Fock term

    dm: density matrix
       <c^\dagger_i c_j>
    
    u_mat: Coulomb tensor
        H_U = (1/2) sum_{ijkl} U_{ijkl} c^\dagger_i c^\dagger_j c_l c_k

    1. Compute asymmetrized Coulomb tensor (asymm_u)
        H_U = (1/4) sum_{ijkl} asymU_{ijkl} c^\dagger_i c^\dagger_j c_l c_k,
        where
           asymU_{ijkl} = -asymU_{ijlk} and
           asymU_{ijkl} = -asymU_{jikl}.
    
    2. Return
        hf_{ij} = sum_{kl} U_{ikjl} <c^\dagger_k c_l>
    """
    nso = dm.shape[0]
    assert dm.shape == (nso, nso)
    assert u_mat.shape == 4*(nso,)

    asymm_u = u_mat.copy()
    asymm_u = 0.5*(asymm_u - asymm_u.transpose((0,1,3,2)))
    asymm_u = 0.5*(asymm_u - asymm_u.transpose((1,0,2,3)))
    asymm_u *= 2
    #for i in range(2):
        #for j in range(2):
            #for k in range(2):
                #for l in range(2):
                    #print(i, j, k, l, asymm_u[i,j,k,l])
    return numpy.einsum('ikjl,kl->ij', asymm_u, dm)


def hf_dc(dm: Dict, u_mat: numpy.ndarray, use_spin_orbit: bool) -> numpy.ndarray:
    """
    Compute HF_DFT or HF_imp from density matrix
    """
    assert isinstance(dm, dict)
    assert isinstance(u_mat, numpy.ndarray)

    if use_spin_orbit:
        assert "ud" in dm.keys()
        return {"ud": hartree_fock_term(dm["ud"], u_mat)}
    else:
        assert "up" in dm.keys() and "down" in dm.keys()
        norb = dm["up"].shape[0]
        nso = 2*norb
        dm_mat = numpy.zeros((2, nso//2, 2, nso//2), dtype=numpy.complex128)
        dm_mat[0,:,0,:] = dm["up"]
        dm_mat[1,:,1,:] = dm["down"]
        dm_mat = dm_mat.reshape(nso, nso)
        hf = hartree_fock_term(dm_mat, u_mat).reshape((2,norb,2,norb))
        return {"up": hf[0,:,0,:], "down": hf[1,:,1,:]}
