import mpi4py
from mpi4py import MPI
import numpy

MPI_TYPE_MAP = {
    'int8': MPI.CHAR,
    'int16': MPI.SHORT,
    'int32': MPI.INT,
    'int64': MPI.LONG,
    'int128': MPI.LONG_LONG,
    'float32': MPI.FLOAT,
    'float64': MPI.DOUBLE,
    'bool': MPI.BOOL,
    'complex128': MPI.DOUBLE_COMPLEX
}

def allgatherv(send_buf):
    """allgatherv for 1D array"""
    MPI = mpi4py.MPI
    comm = MPI.COMM_WORLD
    send_buf = numpy.ascontiguousarray(send_buf).ravel()
    sizes = comm.allgather(send_buf.size)
    offsets = numpy.hstack([0, numpy.cumsum(sizes)])[0:comm.Get_size()]
    recv_buf = numpy.empty(numpy.sum(sizes), dtype=send_buf.dtype)
    MPI_TYPE = MPI_TYPE_MAP[str(send_buf.dtype)]
    comm.Allgatherv([send_buf, MPI_TYPE], [recv_buf, sizes, offsets, MPI_TYPE])
    return recv_buf


def gatherv(send_buf, root=0):
    """gatherv for 1D array"""
    MPI = mpi4py.MPI
    comm = MPI.COMM_WORLD
    send_buf = numpy.ascontiguousarray(send_buf).ravel()
    sizes = comm.allgather(send_buf.size)
    offsets = numpy.hstack([0, numpy.cumsum(sizes)])[0:comm.Get_size()]
    recv_buf = numpy.empty(numpy.sum(sizes), dtype=send_buf.dtype)
    MPI_TYPE = MPI_TYPE_MAP[str(send_buf.dtype)]
    if comm.Get_rank() == root:
        comm.Gatherv([send_buf, MPI_TYPE], [recv_buf, sizes, offsets, MPI_TYPE])
        return recv_buf
    else:
        return None

def split_idx(size, comm_size):
    """Compute sizes and offsets for splitting a 1D array

    Parameters
    ----------
    size : Int
        Length of array
    comm_size : Int
        Number of MPI  processes

    Returns
    -------
    sizes : Int[:]
        Sizes for each MPI processes
    offsets : Int[:]
        Offsets for each MPI processes
    """
    base = size // comm_size
    leftover = int(size % comm_size)

    sizes = numpy.ones(comm_size, dtype=int) * base
    sizes[:leftover] += 1

    offsets = numpy.zeros(comm_size, dtype=int)
    offsets[1:] = numpy.cumsum(sizes)[:-1]

    return sizes, offsets

def get_slice(idx_size, comm=None):
    # get slice for 1D array of size idx_size
    if comm is None:
        comm = mpi4py.MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    sizes, offsets = split_idx(idx_size, size)
    return slice(offsets[rank], offsets[rank]+sizes[rank])


def get_range(idx_size, comm=None):
    # get range for 1D array of size idx_size
    if comm is None:
        comm = mpi4py.MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    sizes, offsets = split_idx(idx_size, size)
    return range(offsets[rank], offsets[rank]+sizes[rank])