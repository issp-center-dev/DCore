from __future__ import print_function

import sys
import os
from pytriqs.applications.dft.converters import *
from pytriqs.archive import *
from pytriqs.utility.h5diff import h5diff
import pytriqs.utility.mpi as mpi

from pytriqs.applications.pydmft.pydmft_pre import pydmft_pre

#
# Execute pydmft_pre.py to generate test.h5
# Then Check the Diff of test.h5 and the reference output (stan_ref.h5))
#
for lattice in ['bethe', 'chain', 'square', 'cubic']:
    input_fname = 'stan_'+lattice+'.in'
    seedname = 'stan_test_' + lattice
    seedname_ref = 'stan_ref_' + lattice

    with open(input_fname, 'w') as f:
        print("[model]", file=f)
        print("t = 1.0", file=f)
        print("U = 4.0", file=f)
        print("lattice = ", lattice, file=f)
        print("seedname = " + seedname, file=f)

    pydmft_pre(input_fname)
    
    print("\n Check Diff of {0}.h5 {1}.h5\n".format(seedname, seedname_ref))
    check = os.system('h5diff -d 1.0e-10 {0}.h5 {1}.h5'.format(seedname, seedname_ref))
    
    if check != 0:
        print("Generated file is different from the refference for lattice = " + lattice)
        sys.exit(-1)
