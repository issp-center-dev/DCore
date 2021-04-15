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

from dcore.numdiff import numdiff
from dcore.openmx2dcore import openmx2dcore

import os

def test_openmx2dcore(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    openmx2dcore("SrVO3", "srvo3")
    numdiff("srvo3_hr.dat", "srvo3_hr_ref.dat")
    numdiff("srvo3_band.dat", "srvo3_band_ref.dat")

    os.chdir(org_dir)