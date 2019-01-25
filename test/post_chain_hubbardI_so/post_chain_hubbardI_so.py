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
from dcore.numdiff import numdiff
from dcore.dcore_post import dcore_post

dcore_post('dmft.ini')

numdiff("test_akw.dat", "test_akw_ref.dat")
numdiff("test_dos.dat", "test_dos_ref.dat")
numdiff("test_momdist.dat", "test_momdist_ref.dat")
numdiff("test_akw0.dat", "test_akw0_ref.dat")
numdiff("test_dos0.dat", "test_dos0_ref.dat")
