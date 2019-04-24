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
from pytriqs.utility.h5diff import h5diff
from dcore.dcore_pre import dcore_pre
#
# Execute dcore_pre.py to generate test.h5
# Then Check the Diff of test.h5 and the reference output (stan_ref.h5))
#
with open("dmft.ini", 'w') as f:
    print("[model]", file=f)
    print("lattice = wannier90", file=f)
    print("seedname = chain", file=f)
    print("nelec = 5.0", file=f)
    print("ncor = 2", file=f)
    print("norb = 3, 4", file=f)
    print("kanamori = [(3.0, 2.0, 0.5), (5.0, 3.0, 1.0)]", file=f)
    print("spin_orbit = True", file=f)
    print("", file=f)
    print("[system]", file=f)
    print("nk0 = 40", file=f)
    print("nk1 = 1", file=f)
    print("nk2 = 1", file=f)

dcore_pre("dmft.ini")

h5diff("chain.h5", "chain_ref.h5")
