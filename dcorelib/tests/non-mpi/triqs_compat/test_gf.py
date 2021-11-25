import pytest
import numpy as np
from dcorelib.sparse_gf.basis import finite_temp_basis
from dcorelib.triqs_compat.gf import *
from dcorelib.triqs_compat.gf.gf import GfImTime
from dcorelib.triqs_compat.gf.tools import *
from dcorelib.triqs_compat.h5 import HDFArchive as HDFArchive2

triqs_available = False
try:
    import triqs.gf as tgf
    from h5 import HDFArchive
    triqs_available = True
except ImportError:
    pass

@pytest.mark.skipif(not triqs_available, reason="TRIQS is not installed.")
def test_gfindices():
    left = ['aaa', 'b']
    right = ['AA', 'B']

    gfindices = tgf.GfIndices([left, right])
    with HDFArchive('gfindices_by_triqs.h5', 'w') as f:
        f['data'] = gfindices

    gfindices_compat = GfIndices([left, right])
    with HDFArchive2('gfindices_by_triqs_compat.h5', 'w') as f:
        f['data'] = gfindices_compat

    for file in ['gfindices_by_triqs.h5', 'gfindices_by_triqs_compat.h5']:
        for HA in [HDFArchive, HDFArchive2]:
            with HA(file, 'r') as f:
                idx = f['data']
                assert list(idx[0]) == left
                assert list(idx[1]) == right


@pytest.mark.skipif(not triqs_available, reason="TRIQS is not installed.")
def test_gf():
    beta = 10.0
    shape = (3,1)
    npoints = 4
    data = np.zeros((2*npoints,) + shape)

    with HDFArchive('gf_triqs.h5', 'w') as f:
        f['gf'] = tgf.GfImFreq(beta=beta, data=data, n_points=npoints)

    with HDFArchive2('gf_triqs_compat.h5', 'w') as f:
        f['gf'] = Gf(beta=beta, data=data)

    # triqs_compat can read the HDF5 file created by TRIQS?
    with HDFArchive2('gf_triqs.h5', 'r') as f:
        _ = f['gf']
    
    # TRIQS can read the HDF5 file created by triqs_compat?
    with HDFArchive('gf_triqs_compat.h5', 'r') as f:
        _ = f['gf/mesh']
        _ = f['gf/indices']
        _ = f['gf']


@pytest.mark.skipif(not triqs_available, reason="TRIQS is not installed.")
def test_block_gf():
    beta = 10.0
    shape = (2,1)
    npoints = 100
    data = np.zeros((2*npoints,) + shape)

    for g in [tgf.GfReFreq(window=(-1,1), n_points=npoints, target_shape=(2,2)),
              tgf.GfImFreq(beta=beta, data=data, n_points=npoints)]:
        tbgf = tgf.BlockGf(name_list=['up', 'dn'], block_list = 2*[g], make_copies = True)

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


def test_tail_fit():
    beta = 10.0
    eps = 1e-7
    basis = finite_temp_basis(beta, 'F', eps=eps)

    h0 = np.array([[1.0, 1j], [-1j, 0.0]])
    nf = h0.shape[0]
    nw = 10000

    giv = GfImFreq(beta=beta, n_points=nw, target_shape=(nf,nf))
    giv << inverse(iOmega_n - h0)

    tail, _ = fit_hermitian_tail(giv, basis)
    
    assert np.abs(h0 - tail[2]).max() < 1e+3*eps

def test_gf_view():
    beta = 10.0
    g = GfImFreq(beta=beta, statistic="F", target_shape=(2,2), n_points=1)
    g1 = GfImFreq(beta=beta, statistic="F", target_shape=(1,1), n_points=1)
    g1.data[...] = 1.0
    g[0,0] << g1
    ref = np.zeros((2,2,2))
    ref[:,0,0] = 1
    assert np.array_equal(g.data,  ref)

def test_block_gf_iter():
    beta = 10.0
    nf = 2
    nw = 10

    bfg = BlockGf(name_list=['up', 'dn'], block_list=
        2*[GfImFreq(beta=beta, n_points=nw, target_shape=(nf,nf))], make_copies=True)
    for ib, (bl, g) in enumerate(bfg):
        g.data[...] = 1.0

    for bl, g in bfg:
        assert (g.data[...] == 1.0).all()
    

def test_delta():
    beta = 10.0
    n_iw = 1000
    D = 2.0 # Half band width
    nso = 1
    #basis = finite_temp_basis(beta, "F", eps=1e-7)
    basis = None

    Delta_iw = GfImFreq(beta=beta, n_points=n_iw, target_shape=(nso, nso))
    Delta_iw << (D/2.0)**2 * SemiCircular(D)

    # Random local transfer matrix
    H0 = np.random.rand(nso, nso)
    H0 = H0 + H0.conjugate().transpose()
    H0[0, 0] = 0.1

    G0_iw = Delta_iw.copy()
    G0_iw << inverse(iOmega_n - H0 - Delta_iw)

    Delta_iw_reconst = delta(G0_iw, basis)

    # Check Delta(iwn)
    diff = Delta_iw_reconst - Delta_iw

    # FIXME: The precision is worse with sparse sampling. Why?
    assert np.all(np.abs(diff.data) < 1e-6)


@pytest.mark.skipif(not triqs_available, reason="TRIQS is not installed.")
def test_conversion():
    beta = 10.0
    n_iw = 10
    nso = 1
    for gclass in [GfImFreq, GfReFreq]:
        gf = gclass(beta=beta, n_points=n_iw, target_shape=(nso, nso))
        gf2 = gf.to_triqs()
        gf3 = Gf.from_triqs(gf2)
        np.testing.assert_array_equal(gf.data, gf3.data)
        np.testing.assert_array_equal(gf.mesh.points, gf3.mesh.points)
