from __future__ import print_function
try:
    from pytriqs import version
except ImportError:
    print('NOT_FOUND', end='')
else:
    print(version.version, end='')