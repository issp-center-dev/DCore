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

import sys
import os
from pytriqs.applications.dft.converters import *
from pytriqs.archive import *
from pytriqs.utility.h5diff import h5diff
import pytriqs.utility.mpi as mpi

from pytriqs.applications.dcore.dcore_pre import dcore_pre

#
# Execute dcore_pre.py to generate test.h5
# Then Check the Diff of test.h5 and the reference output (stan_ref.h5))
#
for lattice in ['bethe', 'chain', 'square', 'cubic']:
    input_fname = 'stan_'+lattice+'.in'
    seedname = 'stan_test_' + lattice
    seedname_ref = 'stan_ref_' + lattice

    with open(input_fname, 'w') as f:
        print("[model]", file=f)
        print("t = 1.0", file=f)
        print("kanamori = [(4.0,0.0,0.0)]", file=f)
        print("lattice = ", lattice, file=f)
        print("seedname = " + seedname, file=f)

    dcore_pre(input_fname)
    
    print("\n Check Diff of {0}.h5 {1}.h5\n".format(seedname, seedname_ref))
    check = os.system('h5diff -d 1.0e-10 {0}.h5 {1}.h5'.format(seedname, seedname_ref))
    
    if check != 0:
        print("Generated file is different from the reference for lattice = " + lattice)
        sys.exit(-1)
