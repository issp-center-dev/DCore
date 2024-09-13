import os, sys
import importlib.util

TRIQS_FOUND = True

triqs_libs = ['triqs', 'triqs_dft_tools']
for l in triqs_libs:
    if not importlib.util.find_spec(l):
        TRIQS_FOUND = False

if "DCORE_TRIQS_COMPAT" in os.environ:
    dtc = os.environ["DCORE_TRIQS_COMPAT"]
    if dtc == "1":
        TRIQS_COMPAT = True
    else:
        TRIQS_COMPAT = False

    if TRIQS_FOUND and TRIQS_COMPAT:
        print("INFO: TRIQS is found but DCORE_TRIQS_COMPAT is set to 1.")
        print("      dcorelib.triqs_compat will be used")

    if not TRIQS_FOUND and not TRIQS_COMPAT:
        print("ERROR: TRIQS is required (DCORE_TRIQS_COMPAT={}) but TRIQS is not found!".format(dtc))
        sys.exit(1)
else:
    TRIQS_COMPAT = not TRIQS_FOUND


if TRIQS_COMPAT:
    from dcorelib.triqs_compat import *
    from dcorelib.triqs_compat import h5
    from dcorelib.triqs_compat.gf import *
    from dcorelib.triqs_compat.gf.gf import make_zero_tail
    from dcorelib.triqs_compat.gf.tools import *
    from dcorelib.triqs_compat.h5 import HDFArchive
    from dcorelib.triqs_compat.utility import *
    from dcorelib.triqs_compat.operators import *
    from dcorelib.triqs_compat import mpi
    from dcorelib.triqs_compat.dft_tools import SumkDFT, SumkDFTTools
    from dcorelib.triqs_compat.plot import mpl_interface
else:
    from triqs.gf.gf import *
    from triqs.gf import *
    from h5 import *
    import h5
    from triqs.utility.h5diff import h5diff, compare, failures
    from triqs.operators.util.op_struct import set_operator_structure
    from triqs.operators.util.U_matrix import *
    from triqs.plot import mpl_interface
    if "OMPI_COMM_WORLD_RANK" in os.environ or "PMI_RANK" in os.environ:
        from triqs.utility import mpi
        from triqs_dft_tools import SumkDFT, SumkDFTTools
    else:
        from .backend import _triqs_mpi as mpi

    print("INFO: TRIQS library is used")