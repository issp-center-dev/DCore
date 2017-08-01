from __future__ import print_function

import sys
import os
from pytriqs.applications.dft.converters import *
from pytriqs.archive import *
from pytriqs.utility.h5diff import h5diff
import pytriqs.utility.mpi as mpi

from pytriqs.applications.pydmft.standard import standard

#
# Execute standard.py to generate test.h5
# Then Check the Diff of test.h5 and the reference output (stan_ref.h5))
#
standard('stan.in')

print("\n Check Diff of test.h5 stan_ref.h5\n")
check = os.system('h5diff test.h5 stan_ref.h5')

if check == 0:
    print("\n Generated file is the same as the refference.\n")
    sys.exit(0)
else:
    print("\n Generated file is different from the refference.\n")
    sys.exit(-1)
