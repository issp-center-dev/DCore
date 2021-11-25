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
from h5 import *
from dcorelib.triqs_compat.gf import *
from dcorelib.triqs_compat.dft_tools.sumk_dft import *
import dcorelib.triqs_compat.mpi
from dcorelib.triqs_compat.operators.util import set_operator_structure
from dcorelib.triqs_compat.utility.comparison_tests import *
from dcorelib.triqs_compat.utility.h5diff import h5diff

def test_Gloc(request):
    os.chdir(request.fspath.dirname)
    # Basic input parameters
    beta = 40

    # Init the SumK class
    SK=SumkDFT(hdf_file='SrVO3.ref.h5',use_dft_blocks=True)

    num_orbitals = SK.corr_shells[0]['dim']
    l = SK.corr_shells[0]['l']
    spin_names = ['down','up']
    orb_names = ['%s'%i for i in range(num_orbitals)]
    orb_hybridized = False

    gf_struct = set_operator_structure(spin_names,orb_names,orb_hybridized)
    glist = [ GfImFreq(indices=inner,beta=beta) for block,inner in gf_struct]
    Sigma_iw = BlockGf(name_list = [block for block,inner in gf_struct], block_list = glist, make_copies = False)

    SK.set_Sigma([Sigma_iw])
    Gloc = SK.extract_G_loc()

    if mpi.is_master_node():
        with HDFArchive('srvo3_Gloc.out.h5','w') as ar:
            ar['Gloc'] = Gloc[0]

        h5diff("srvo3_Gloc.out.h5","srvo3_Gloc.ref.h5", precision=1e-14)
