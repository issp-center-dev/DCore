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

from triqs.utility.h5diff import h5diff
from dcore.dcore_pre import dcore_pre
import os

#
# Execute dcore_pre.py to generate test.h5
# Then Check the Diff of test.h5 and the reference output (stan_ref.h5))
#
def test_pre_wannier(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    with open("nis.ini", 'w') as f:
        print("[model]", file=f)
        print("lattice = wannier90", file=f)
        print("seedname = nis", file=f)
        print("nelec = 24.0", file=f)
        print("ncor = 2", file=f)
        print("norb = [5, 5]", file=f)
        print("interaction = slater_uj", file=f)
        print("slater_uj = [(2, 1.0, 0.0), (2, 1.0, 0.0)]", file=f)
        print("nk0 = 4", file=f)
        print("nk1 = 4", file=f)
        print("nk2 = 3", file=f)
        print("", file=f)
    
    dcore_pre("nis.ini")
    
    h5diff("nis.h5", "nis_ref.h5")

    os.chdir(org_dir)