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

from dcore.tools import spin_moments_sh
import numpy
import scipy

def test_spin_moments_sh():
    # (Sx, Sy, Sz) = (0, 0, 1/2)
    dm_corr_sh = [{'ud': numpy.array([[1, 0], [0, 0]])}]
    assert numpy.allclose(spin_moments_sh(dm_corr_sh)[0], numpy.array([0.0, 0.0, 0.5]) )
    
    dm_corr_sh = [{'up': numpy.array([[1]]), 'down': numpy.array([[0]])}]
    assert numpy.allclose(spin_moments_sh(dm_corr_sh)[0], numpy.array([0.0, 0.0, 0.5]) )
    
    # (Sx, Sy, Sz) = (1/2, 0, 0)
    dm_corr_sh = [{'ud': 0.5*numpy.ones((2,2))}]
    assert numpy.allclose(spin_moments_sh(dm_corr_sh)[0], numpy.array([0.5, 0.0, 0.0]) )

    # Two-orbital case
    norb = 2
    dm_mat = numpy.zeros((2,norb,2,norb))
    for iorb in range(norb):
        dm_mat[:,iorb,:,iorb] = numpy.array([[1, 0], [0, 0]])
    dm_corr_sh = [{'ud': dm_mat.reshape(2*norb,2*norb)}]
    assert numpy.allclose(spin_moments_sh(dm_corr_sh)[0], numpy.array([0, 0, norb*0.5]) )

test_spin_moments_sh()
