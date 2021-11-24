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

from dcore.dcore_pre import dcore_pre
from dcore.backend import h5diff
import shutil, os

#
# Execute dcore_pre.py to generate test.h5
# Then Check the Diff of test.h5 and the reference output (stan_ref.h5))
#

def test_pre_external(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    input_ini = "external.in"
    shutil.copy("external.h5.in", "external.h5")
    
    with open(input_ini, 'w') as f:
        print("[model]", file=f)
        print("lattice = external", file=f)
        print("seedname = external", file=f)
        print("kanamori = [(4.0,0.0,0.0)]", file=f)
    
    dcore_pre(input_ini)
    
    h5diff("external.h5", "external_ref.h5")

    os.chdir(org_dir)
