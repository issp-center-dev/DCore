import numpy as np
from dcore.backend.sparse_gf.transform import fourier

from dcore.backend.triqs_compat.gf import *
from dcore.backend.triqs_compat.gf.gf import GfImTime
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

def test_transform():
    beta = 10.0
    h0 = np.array([[1.0, 1j], [-1j, 0.0]])
    nf = h0.shape[0]
    assert np.abs(h0-h0.T.conj()).max() < 1e-10
    nw = 1000
    eps = 1e-5

    # Compute Matsubara Green's function
    iv = 1J*(2*np.arange(-nw,nw)+1)*np.pi/beta
    iv_mat = np.einsum('w,ij->wij', iv, np.identity(nf))
    data = np.linalg.inv(iv_mat - h0[None,:,:])
    giv = GfImFreq(data=data, beta=beta)

    # Create intermediate object `ft` from Matsubara
    ft = fourier(giv)

    # Reconstruct giv
    giv_reconst = giv.copy()
    giv_reconst << ft
    assert np.abs(giv.data-giv_reconst.data).max() < eps

    # Comppute imaginary-time Green's function
    ntau = 10000
    taus = np.linspace(0, beta, ntau)
    evals, evecs = np.linalg.eigh(h0)
    gtau_diag_data = np.zeros((ntau, nf, nf), dtype=np.complex128)
    for ie in range(nf):
        gtau_diag_data[:,ie,ie] = -np.exp(-taus*evals[ie])/(1+np.exp(-beta*evals[ie]))
    gtau_diag = GfImTime(data=gtau_diag_data, beta=beta)
    gtau = gtau_diag.copy()
    gtau.from_L_G_R(evecs, gtau_diag, evecs.T.conj())
    
    # Create intermediate object `ft` from tau
    ft2 = fourier(gtau)
    giv_reconst2 = giv.copy()
    giv_reconst2 << ft2
    assert np.abs(giv.data-giv_reconst2.data).max() < eps

def test_lazy():
    beta = 10.0
    eps = 1e-5

    h0 = np.array([[1.0, 1j], [-1j, 0.0]])
    nf = h0.shape[0]
    assert np.abs(h0-h0.T.conj()).max() < 1e-10
    nw = 1000

    # Compute Matsubara Green's function
    iv = 1J*(2*np.arange(-nw,nw)+1)*np.pi/beta
    iv_mat = np.einsum('w,ij->wij', iv, np.identity(nf))
    data = np.linalg.inv(iv_mat - h0[None,:,:])
    giv_ref = GfImFreq(data=data, beta=beta)

    giv_from_lazy = giv_ref.copy()
    giv_from_lazy.zero()
    giv_from_lazy << inverse(iOmega_n - h0)
    assert np.abs(giv_ref.data-giv_from_lazy.data).max() < eps

    h0_ = giv_ref.copy()
    h0_ << h0
    giv_from_lazy << inverse(iOmega_n - h0_)
    assert np.abs(giv_ref.data-giv_from_lazy.data).max() < eps