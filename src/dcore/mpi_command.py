import os

default_mpi_command = "mpirun -np #"

if 'DCORE_MPIRUN_COMMAND' in os.environ:
    default_mpi_command = os.environ['DCORE_MPIRUN_COMMAND']
