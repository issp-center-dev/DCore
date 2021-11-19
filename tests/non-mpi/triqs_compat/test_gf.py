import numpy as np

from dcore.backend.triqs_compat.gf import *
from dcore.backend.triqs_compat.h5 import HDFArchive as HDFArchive2

import triqs.gf as tgf
from h5 import HDFArchive

def test_gfindices():
    left = ['aaa', 'b']
    right = ['AA', 'B']

    gfindices = tgf.GfIndices([left, right])
    with HDFArchive('gfindices_by_triqs.h5', 'w') as f:
        f['data'] = gfindices

    gfindices_compat = GfIndices([left, right])
    with HDFArchive2('gfindices_by_triqs_compat.h5', 'w') as f:
        f['data'] = gfindices_compat

    # TRIQS can read the HDF5 file created by triqs_compat?
    with HDFArchive('gfindices_by_triqs_compat.h5', 'r') as f:
        idx = f['data']
        assert list(idx[0]) == left
        assert list(idx[1]) == right

    # triqs_compat can read the HDF5 file created by TRIQS?
    with HDFArchive2('gfindices_by_triqs.h5', 'r') as f:
        idx = f['data']
        assert list(idx[0]) == left
        assert list(idx[1]) == right

def test_gf():
    beta = 10.0
    shape = (3,1)
    npoints = 4
    data = np.zeros((2*npoints,) + shape)

    gf = tgf.GfImFreq(beta=beta, data=data, n_points=npoints)
    with HDFArchive('gf.h5', 'w') as f:
        f['gf'] = gf

    gf2 = Gf(beta=beta, data=data)
    with HDFArchive2('gf2.h5', 'w') as f:
        f['gf'] = gf2

    # TRIQS can read the HDF5 file created by triqs_compat?
    with HDFArchive('gf2.h5', 'r') as f:
        gf = f['gf']

    # triqs_compat can read the HDF5 file created by TRIQS?
    with HDFArchive2('gf.h5', 'r') as f:
        gf = f['gf']
    

def test_block_gf():
    beta = 10.0
    shape = (2,1)
    npoints = 100
    data = np.zeros((2*npoints,) + shape)

    tbgf = tgf.BlockGf(name_list=['up', 'dn'],
        block_list = 2*[tgf.GfImFreq(beta=beta, data=data, n_points=npoints)],
        make_copies = True)
    with HDFArchive('bgf_by_triqs.h5', 'w') as f:
        f['data'] = tbgf

    bgf = BlockGf(name_list=['up', 'dn'],
        block_list = 2*[Gf(beta=beta, data=data)],
        make_copies = True)

    with HDFArchive2('bgf_by_compat.h5', 'w') as f:
        f['data'] = bgf

    with HDFArchive('bgf_by_compat.h5', 'r') as f:
        f['data']

    with HDFArchive2('bgf_by_triqs.h5', 'r') as f:
        f['data']
        print(f['data'])

    