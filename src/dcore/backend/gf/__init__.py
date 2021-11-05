from .. import _backend_module

if _backend_module == 'TRIQS_COMPAT':
    from ..triqs_compat.gf import *
elif _backend_module == 'TRIQS':
    from triqs.gf import *