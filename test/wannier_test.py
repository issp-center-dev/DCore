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
with open("nis.ini", 'w') as f:
    print("[model]", file=f)
    print("U = 1.0", file=f)
    print("lattice = wannier90", file=f)
    print("seedname = nis", file=f)
    print("nelec = 24.0", file=f)
    print("ncor = 2", file=f)
    print("cshell = [(2, 5), (2, 5)]", file=f)
    print("", file=f)
    print("[system]", file=f)
    print("nk0 = 4", file=f)
    print("nk1 = 4", file=f)
    print("nk2 = 3", file=f)

pydmft_pre("nis.ini")
    
print("\n Check Diff of nis.h5 nis_ref.h5\n")
check = os.system('h5diff -d 1.0e-10 nis.h5 nis_ref.h5')
    
if check != 0:
    print("Generated file and reference file are different.")
    sys.exit(-1)
