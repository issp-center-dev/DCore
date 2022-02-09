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


import glob
import re
import os
from dcore.tools import h5diff
from dcore.numdiff import numdiff
from dcore.dcore_mpicheck import dcore_mpicheck
from dcore.dcore_pre import dcore_pre
from dcore.dcore import dcore
from dcore.dcore_check import dcore_check
from dcore.dcore_post import dcore_post
from dcore.dcore_pade import dcore_pade


def test_chain_hubbardI(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    seedname = 'test'
    dcore_mpicheck('dmft.ini')
    dcore_pre('dmft.ini')
    dcore('dmft.ini')
    dcore_check('dmft.ini', './check', 'eps', 10000)
    dcore_pade(seedname)
    dcore_post('dmft.ini')
    
    data_files = glob.glob('./ref/*')
    
    for path in data_files:
        base_name = os.path.basename(path)
        print("base_nam,e ", base_name)
        if base_name == seedname + '.h5':
            h5diff(base_name, path)
        elif base_name == seedname + '.out.h5':
            h5diff(base_name, path, "dmft_out/Sigma_iw", precision=1e-2)
        elif not re.search('.dat$', base_name) is None:
            numdiff(base_name, path, 1e-2)
        else:
            raise RuntimeError("Uknown how to check " + base_name)

    os.chdir(org_dir)
