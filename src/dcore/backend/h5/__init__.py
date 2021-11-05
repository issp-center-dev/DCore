from .. import _backend_module

if _backend_module == 'TRIQS_COMPAT':
    from ..triqs_compat.h5 import *
elif _backend_module == 'TRIQS':
    from h5 import *