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
#

import os

def test_pre_respack_so(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    dcore_pre("dmft.ini")
    h5diff("test.h5", "test_ref.h5")

    os.chdir(org_dir)