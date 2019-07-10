from pytriqs import version

triqs_major_version = int(version.version.split('.')[0])

if triqs_major_version == 1:
    from pytriqs.gf.local import *
elif triqs_major_version == 2:
    from pytriqs.gf import *
else:
    raise RuntimeError("Unknown triqs major version!")
