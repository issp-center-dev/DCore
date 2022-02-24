from dcore.mpi import bcast
from dcore._dispatcher import mpi
import numpy

def test_bcast_split_transfer():
    max_bytes = 100000000
    for size in [int(0.1*max_bytes)//8, max_bytes//8]:
        send_buff = numpy.arange(size) if mpi.is_master_node() else None
        recv_buff = bcast(send_buff, max_bytes=max_bytes)
        numpy.testing.assert_array_equal(numpy.arange(size), recv_buff)