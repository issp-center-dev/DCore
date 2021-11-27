
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


from dcorelib.triqs_compat.dft_tools.converters import *
from dcorelib.triqs_compat.h5 import *
from dcorelib.triqs_compat.utility.h5diff import h5diff

def test_w90_convert():
    Converter = Wannier90Converter(seedname='LaVO3-Pnma',hdf_filename='w90_convert.out.h5')
    Converter.convert_dft_input()
    h5diff("w90_convert.out.h5","w90_convert.ref.h5") 
