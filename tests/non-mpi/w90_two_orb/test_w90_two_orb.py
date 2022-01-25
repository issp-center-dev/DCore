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

import os
import numpy
from dcore.dcore_pre import dcore_pre
from dcore._dispatcher import HDFArchive


def test_dcore_pre(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    dcore_pre('dmft_square.ini')
    
    nk1 = 32
    nk2 = 32
    nk3 = 1
    norb = 2
    onsite_ene = [-0.2, 0.2] # CF slit
    bw = [1, 0.5] # Band width

    #0  0  1 1  -0.2 0.0
    #0  0 0  1 2   0.0 0.0
    #0  0 0  2 1   0.0 0.0
    #0  0 0  2 2   0.2 0.0

    with HDFArchive("square_2orb.h5", "r") as f:
        hk = f["dft_input"]["hopping"].reshape((nk1, nk2, nk3, norb, norb))

    for ik1 in range(nk1):
        for ik2 in range(nk2):
            k1 = ik1 * 2*numpy.pi/nk1
            k2 = ik2 * 2*numpy.pi/nk2
            for iorb in range(2):
                hk_ref = bw[iorb] * (- 2*numpy.cos(k1) - 2*numpy.cos(k2)) + onsite_ene[iorb]
                assert numpy.abs(hk_ref-hk[ik1, ik2, 0, iorb, iorb]) < 1e-10