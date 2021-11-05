import os

from .triqs_compat import *  # REMOVE THIS AT SOME POINT

_backend_module = 'TRIQS_COMPAT'

#if 'DCORE_TRIQS_COMPAT' in os.environ and os.environ['DCORE_TRIQS_COMPAT'] == 1:
    #from .triqs_compat import *
#else:
    #from .triqs import *