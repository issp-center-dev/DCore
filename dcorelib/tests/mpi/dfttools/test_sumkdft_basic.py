
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

import os
from dcorelib.triqs_compat.h5 import *
from dcorelib.triqs_compat.dft_tools.sumk_dft_tools import SumkDFTTools
from dcorelib.triqs_compat import mpi
from dcorelib.triqs_compat.utility.comparison_tests import *
from dcorelib.triqs_compat.utility.h5diff import h5diff


def test_svo(request):
    os.chdir(request.fspath.dirname)
    SK = SumkDFTTools(hdf_file = 'SrVO3.ref.h5')
    dm = SK.density_matrix(method = 'using_gf', beta = 40)
    dm_pc = SK.partial_charges(beta=40,with_Sigma=False,with_dc=False)
    if mpi.is_master_node():
        with HDFArchive('sumkdft_basic.out.h5','w') as ar:
            ar['dm'] = dm
            ar['dm_pc'] = dm_pc
        h5diff('sumkdft_basic.out.h5','sumkdft_basic.ref.h5', precision=1e-5)
