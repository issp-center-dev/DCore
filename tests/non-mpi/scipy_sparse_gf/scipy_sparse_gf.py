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
# Load the standalone solver without importing the backend-dependent solver registry.
import importlib.util
from pathlib import Path
import sys
from unittest.mock import patch

_solver_dir = Path(__file__).resolve().parents[3] / 'src/dcore/impurity_solvers'


def _load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, _solver_dir / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_lanczos = _load_module('lanczos', 'lanczos.py')
with patch.dict(sys.modules, {'dcore.impurity_solvers.lanczos': _lanczos}):
    _solver = _load_module('scipy_sparse_main', 'scipy_sparse_main.py')
calc_gf_Lehmann = _solver.calc_gf_Lehmann
_dense_eigh = _solver._dense_eigh


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


@pytest.mark.parametrize("driver", ["evr", "evd", "ev", "evx"])
def test_dense_eigh_cpu_matches_scipy(driver):
    import scipy.linalg, scipy.sparse as sp
    rng = np.random.default_rng(0)
    A = rng.standard_normal((8, 8)) + 1j*rng.standard_normal((8, 8))
    H = sp.csr_matrix(A + A.conj().T)          # Hermitian
    w, v = _dense_eigh(H, np, driver=driver)
    w_ref = scipy.linalg.eigh(H.toarray(), eigvals_only=True)
    assert np.allclose(np.sort(w), np.sort(w_ref), atol=1e-10)
    # eigenpairs reconstruct H
    assert np.allclose((v * w) @ v.conj().T, H.toarray(), atol=1e-9)


@pytest.mark.parametrize("driver", ["evr", "evd"])
def test_dense_eigh_degenerate_spectrum(driver):
    # Highly degenerate spectrum, as in particle-hole symmetric ED blocks
    import scipy.sparse as sp
    rng = np.random.default_rng(3)
    Q, _ = np.linalg.qr(rng.standard_normal((40, 40)))
    w_ref = np.repeat([-1.0, 0.0, 2.0, 5.0], 10)
    H = sp.csr_matrix((Q * w_ref) @ Q.T)
    w, v = _dense_eigh(H, np, driver=driver)
    np.testing.assert_allclose(w, w_ref, atol=1e-10)
    np.testing.assert_allclose(v.T @ v, np.eye(40), atol=1e-10)
    np.testing.assert_allclose((v * w) @ v.T, H.toarray(), atol=1e-10)


def test_invalid_eigh_driver_rejected(tmp_path, monkeypatch):
    import json
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['scipy_sparse_main.py', 'input.json'])
    Path('input.json').write_text(json.dumps(dict(
        n_flavors=2, n_sites=1, beta=5.0, n_eigen=4, n_iw=13, flag_spin_conserve=1,
        dim_full_diag=10, particle_numbers='all', weight_threshold=0.0, ncv=None,
        eigen_solver='eigsh', gf_solver='bicgstab', check_n_eigen=True,
        check_orthonormality=True, file_h0='h0.npy', file_umat='umat.npy',
        eigh_driver='gvd')))
    with pytest.raises(ValueError, match="eigh_driver"):
        _solver.main()


def test_dense_eigh_default_driver_without_driver_argument(monkeypatch):
    # SciPy < 1.5: scipy.linalg.eigh has no driver argument
    import scipy.linalg, scipy.sparse as sp
    orig = scipy.linalg.eigh
    monkeypatch.setattr(scipy.linalg, 'eigh', lambda a: orig(a))
    H = sp.csr_matrix(np.diag([3.0, 1.0, 2.0]))
    w, _ = _dense_eigh(H, np)
    np.testing.assert_allclose(w, [1.0, 2.0, 3.0])


def test_eigh_driver_rejected_without_driver_argument(tmp_path, monkeypatch):
    import json, scipy.linalg
    orig = scipy.linalg.eigh
    monkeypatch.setattr(scipy.linalg, 'eigh', lambda a: orig(a))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['scipy_sparse_main.py', 'input.json'])
    Path('input.json').write_text(json.dumps(dict(
        n_flavors=2, n_sites=1, beta=5.0, n_eigen=4, n_iw=13, flag_spin_conserve=1,
        dim_full_diag=10, particle_numbers='all', weight_threshold=0.0, ncv=None,
        eigen_solver='eigsh', gf_solver='bicgstab', check_n_eigen=True,
        check_orthonormality=True, file_h0='h0.npy', file_umat='umat.npy',
        eigh_driver='evd')))
    with pytest.raises(ValueError, match="SciPy >= 1.5"):
        _solver.main()


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


@pytest.mark.parametrize("pm", [+1, -1])
@pytest.mark.parametrize("spin_conserve", [True, False])
@pytest.mark.parametrize("device", [False, True])
def test_lehmann_bounds_frequency_workspace(monkeypatch, pm, spin_conserve, device):
    iws, Cdag, eigvec, evx, vvx = _rand_case(seed=19, niw=13)
    monkeypatch.setattr(_solver, '_LEHMANN_CHUNK_BYTES', 5 * evx.size * 16, raising=False)
    shapes = []
    original = np.einsum

    def einsum(expression, ci, ene, cj):
        shapes.append(ene.shape)
        return original(expression, ci, ene, cj)

    ref = _reference_lehmann(iws, Cdag, spin_conserve, eigvec, 0.3, evx, vvx, pm)
    monkeypatch.setattr(np, 'einsum', einsum)
    got = calc_gf_Lehmann(iws, Cdag, spin_conserve, eigvec, 0.3, evx, vvx, pm,
                         xp=_FakeCupy() if device else np)
    np.testing.assert_allclose(got, ref, atol=1e-12)
    assert shapes == [(5, 6), (5, 6), (3, 6)]


@pytest.mark.parametrize('device', [False, True])
def test_atomic_green_function(tmp_path, monkeypatch, device):
    import json
    from dcore import gpu
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['scipy_sparse_main.py', 'input.json'])
    requested = []

    def backend(use_gpu):
        requested.append(use_gpu)
        return (_FakeCupy(), True) if use_gpu else (np, False)

    monkeypatch.setattr(gpu, 'get_backend', backend)
    np.save('h0.npy', -np.eye(2))
    umat = np.zeros((2, 2, 2, 2))
    umat[0, 1, 0, 1] = umat[1, 0, 1, 0] = 2.0
    np.save('umat.npy', umat)
    params = dict(n_flavors=2, n_sites=1, beta=5.0, n_eigen=4, n_iw=13,
                  flag_spin_conserve=1, dim_full_diag=10, particle_numbers='all',
                  weight_threshold=0.0, ncv=None, eigen_solver='eigsh',
                  gf_solver='bicgstab', check_n_eigen=True, check_orthonormality=True,
                  file_h0='h0.npy', file_umat='umat.npy', gpu=device)
    Path('input.json').write_text(json.dumps(params))
    _solver.main()
    iw = 1j * (2 * np.arange(13) + 1) * np.pi / params['beta']
    expected = 0.5 / (iw - 1) + 0.5 / (iw + 1)
    gf = np.load('gf.npy')
    np.testing.assert_allclose(gf[0, 0], expected, atol=1e-12)
    np.testing.assert_allclose(gf[1, 1], expected, atol=1e-12)
    np.testing.assert_allclose(gf[0, 1], 0, atol=1e-12)
    np.testing.assert_allclose(gf[1, 0], 0, atol=1e-12)
    assert requested == [device]
