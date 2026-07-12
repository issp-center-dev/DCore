#
# DCore -- Integrated DMFT software for correlated electrons
# Copyright (C) 2017 The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import numpy as np
import scipy.sparse as sp
import pytest
from dcore.impurity_solvers.scipy_sparse_main import calc_gf_Lehmann


def _reference_lehmann(iws, Cdag, spin_conserve, eigvec, E_n, evx, vvx, pm):
    nf = Cdag.size; niw = iws.size; dim = evx.size
    gf = np.zeros((nf, nf, niw), dtype=complex)
    cdag_im = np.empty((nf, dim), dtype=complex)
    for i in range(nf):
        cdag_im[i] = vvx.conj().T @ Cdag[i] @ eigvec
    for l, iw in enumerate(iws):
        den = 1/(iw - evx + E_n) if pm == +1 else 1/(iw + evx - E_n)
        for i, j in np.ndindex(nf, nf):
            if spin_conserve and (i // (nf//2) != j // (nf//2)):
                continue
            gf[i, j, l] = np.einsum("m,m,m", cdag_im[i].conj(), den, cdag_im[j])
    return gf


def _rand_case(seed, nf=4, dim=6, niw=5):
    rng = np.random.default_rng(seed)
    def cx(shape):
        return rng.standard_normal(shape) + 1j*rng.standard_normal(shape)
    Cdag = np.empty(nf, dtype=object)
    for i in range(nf):
        Cdag[i] = sp.csr_matrix(cx((dim, dim)))
    eigvec = cx(dim)
    evx = rng.standard_normal(dim)
    vvx = cx((dim, dim))
    iws = 1j * (2*np.arange(niw)+1) * np.pi / 10.0
    return iws, Cdag, eigvec, evx, vvx


@pytest.mark.parametrize("pm", [+1, -1])
@pytest.mark.parametrize("spin_conserve", [True, False])
def test_lehmann_matches_reference(pm, spin_conserve):
    iws, Cdag, eigvec, evx, vvx = _rand_case(
        seed=(0 if pm == +1 else 1) + (10 if spin_conserve else 0)
    )
    ref = _reference_lehmann(iws, Cdag, spin_conserve, eigvec, 0.3, evx, vvx, pm)
    got = calc_gf_Lehmann(iws, Cdag, spin_conserve, eigvec, 0.3, evx, vvx, pm)
    assert got.shape == ref.shape
    assert np.allclose(got, ref, atol=1e-12)


def test_dense_eigh_cpu_matches_scipy():
    import scipy.linalg, scipy.sparse as sp
    from dcore.impurity_solvers.scipy_sparse_main import _dense_eigh
    rng = np.random.default_rng(0)
    A = rng.standard_normal((8, 8)) + 1j*rng.standard_normal((8, 8))
    H = sp.csr_matrix(A + A.conj().T)          # Hermitian
    w, v = _dense_eigh(H, np)
    w_ref = scipy.linalg.eigh(H.toarray(), eigvals_only=True)
    assert np.allclose(np.sort(w), np.sort(w_ref), atol=1e-10)
    # eigenpairs reconstruct H
    assert np.allclose((v * w) @ v.conj().T, H.toarray(), atol=1e-9)


class _FakeCupy:
    # minimal shim: behaves like numpy but is a distinct module identity,
    # and provides asnumpy, so the xp-is-not-np branches execute.
    def __getattr__(self, k):
        import numpy as _np
        return getattr(_np, k)
    @staticmethod
    def asnumpy(a):
        import numpy as _np
        return _np.asarray(a)


def test_lehmann_xp_branch_executes():
    iws, Cdag, eigvec, evx, vvx = _rand_case(seed=7)
    ref = _reference_lehmann(iws, Cdag, False, eigvec, 0.1, evx, vvx, +1)
    got = calc_gf_Lehmann(iws, Cdag, False, eigvec, 0.1, evx, vvx, +1, xp=_FakeCupy())
    assert np.allclose(got, ref, atol=1e-12)


def test_gpu_matches_cpu_if_available():
    cupy = pytest.importorskip("cupy")
    try:
        if cupy.cuda.runtime.getDeviceCount() < 1:
            pytest.skip("no CUDA device")
    except Exception:
        pytest.skip("no usable CUDA device")
    from dcore.impurity_solvers.scipy_sparse_main import _dense_eigh
    import scipy.sparse as sp
    rng = np.random.default_rng(1)
    A = rng.standard_normal((32, 32)) + 1j*rng.standard_normal((32, 32))
    H = sp.csr_matrix(A + A.conj().T)
    w_cpu, _ = _dense_eigh(H, np)
    w_gpu, _ = _dense_eigh(H, cupy)
    assert np.allclose(np.sort(w_cpu), np.sort(w_gpu), rtol=1e-10, atol=1e-10)
    iws, Cdag, eigvec, evx, vvx = _rand_case(seed=2)
    g_cpu = calc_gf_Lehmann(iws, Cdag, False, eigvec, 0.2, evx, vvx, +1, xp=np)
    g_gpu = calc_gf_Lehmann(iws, Cdag, False, eigvec, 0.2, evx, vvx, +1, xp=cupy)
    assert np.allclose(g_cpu, g_gpu, rtol=1e-10, atol=1e-12)
