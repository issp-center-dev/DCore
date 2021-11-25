import os, sys

if 'triqs.utility.mpi' in sys.modules:
    raise RuntimeError("Error: MPI must not be imported in a non-MPI module! This indicates a bug in DCore.")

if 'DCORE_TRIQS_COMPAT' in os.environ and int(os.environ['DCORE_TRIQS_COMPAT']) == 1:
    from dcore.backend.triqs_compat import *
    from dcore.backend.triqs_compat import h5
    from dcore.backend.triqs_compat.gf import *
    from dcore.backend.triqs_compat.gf.tools import *
    from dcore.backend.triqs_compat.h5 import HDFArchive
    from dcore.backend.triqs_compat.utility import *
    from dcore.backend.triqs_compat.operators import *
    from dcore.backend.triqs_compat import mpi
    from dcore.backend.triqs_compat.dft_tools import SumkDFT, SumkDFTTools
else:
    from triqs.gf.gf import *
    from triqs.gf import *
    from h5 import *
    import h5
    from triqs.utility.h5diff import h5diff, compare, failures
    from triqs.operators.util.op_struct import set_operator_structure
    from triqs.operators.util.U_matrix import *
    if "OMPI_COMM_WORLD_RANK" in os.environ or "PMI_RANK" in os.environ:
        from triqs.utility import mpi
        from triqs_dft_tools import SumkDFT, SumkDFTTools