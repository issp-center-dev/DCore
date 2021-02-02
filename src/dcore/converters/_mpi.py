import sys

def is_mpi_loaded():
    return 'triqs.utility.mpi' in sys.modules

def is_master_node():
    if is_mpi_loaded():
        import triqs.utility.mpi as mpi
        return mpi.is_master_node()
    else:
        return True

def report(message):
    if is_mpi_loaded():
        import triqs.utility.mpi as mpi
        return mpi.report(message)
    else:
        print(message)