import os, sys
import importlib.util

if "DCORE_TRIQS_COMPAT" in os.environ and os.environ["DCORE_TRIQS_COMPAT"] == "0":
    TRIQS_COMPAT = False
else:
    TRIQS_COMPAT = True

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
    triqs_libs = ['triqs', 'triqs_dft_tools']
    for l in triqs_libs:
        if not importlib.util.find_spec(l):
            sys.exit("ERROR: TRIQS library is not found. Please install TRIQS library or set DCORE_TRIQS_COMPAT=1")

    from triqs.version import version as triqs_version

    if not "2" < triqs_version[0] < "4":
        sys.exit("ERROR: TRIQS version {} is not supported. Please install TRIQS version 3.x".format(triqs_version))
    
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

    print("INFO: TRIQS library {} is used".format(triqs_version))