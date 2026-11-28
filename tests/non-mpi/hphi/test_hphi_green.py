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
"""
Reconstruction of G_ij from the spectra of excitation operators (HPhi solver).
Runs without HPhi: the spectra are built from known particle/hole parts.
"""
import itertools
import numpy as np

from dcore.impurity_solvers.hphi_spectrum import calc_one_body_green


def _spectra(P, H):
    """
    Spectra G_X = <X^+ R X> as computed by HPhi (coefficients applied to the ket):
      ex = 1: X = c^+_a + lam c^+_b  ->  P_aa + |lam|^2 P_bb + lam P_ab + conj(lam) P_ba
      ex = 0: X = c_a + lam c_b      ->  H_aa + |lam|^2 H_bb + lam H_ba + conj(lam) H_ab
    with P_ab = <c_a R c^+_b>, H_ab = <c^+_b R c_a>; lam = 1j (flg = 0) or 1 (flg = 1).
    """
    n_site, n_sigma = P.shape[0], P.shape[1]
    shape = P.shape[:4] + (2, 2) + P.shape[4:]
    core = np.zeros(shape, dtype=complex)
    for (i, si), (j, sj) in itertools.product(itertools.product(range(n_site), range(n_sigma)), repeat=2):
        for flg, lam in [(0, 1j), (1, 1.0)]:
            if (i, si) == (j, sj):
                core[i, si, j, sj, flg, 1] = P[i, si, i, si]
                core[i, si, j, sj, flg, 0] = H[i, si, i, si]
                continue
            core[i, si, j, sj, flg, 1] = (P[i, si, i, si] + abs(lam) ** 2 * P[j, sj, j, sj]
                                          + lam * P[i, si, j, sj] + np.conj(lam) * P[j, sj, i, si])
            core[i, si, j, sj, flg, 0] = (H[i, si, i, si] + abs(lam) ** 2 * H[j, sj, j, sj]
                                          + lam * H[j, sj, i, si] + np.conj(lam) * H[i, si, j, sj])
    return core


def test_reconstruction_random():
    rng = np.random.default_rng(0)
    shape = (2, 2, 2, 2, 1, 5)   # (n_site, n_sigma, n_site, n_sigma, n_T, n_omega)
    P = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    H = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    G = calc_one_body_green(_spectra(P, H))
    assert np.allclose(G, P + H, atol=1e-13)


def test_reconstruction_lehmann():
    """Physical G of a random 2-site fermion model (4 spin-orbitals), Lehmann representation."""
    rng = np.random.default_rng(1)
    n_so = 4
    dim = 2 ** n_so
    cd = np.array([[0, 0], [1, 0]], dtype=float)
    ops = []
    for k in range(n_so):
        m = np.ones((1, 1))
        for f in [np.eye(2)] * k + [cd] + [np.diag([1.0, -1.0])] * (n_so - k - 1):
            m = np.kron(m, f)
        ops.append(m)
    A = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    Ham = A + A.conj().T
    E, V = np.linalg.eigh(Ham)
    beta = 2.0
    w = np.exp(-beta * (E - E[0]))
    w /= w.sum()
    iw = 1j * (2 * np.arange(4) + 1) * np.pi / beta
    cdag = [V.conj().T @ o @ V for o in ops]
    c = [m.conj().T for m in cdag]
    # flavor (site, sigma) -> index 2*site + sigma
    P = np.zeros((2, 2, 2, 2, 1, iw.size), dtype=complex)
    H = np.zeros_like(P)
    for (i, si), (j, sj) in itertools.product(itertools.product(range(2), range(2)), repeat=2):
        a, b = 2 * i + si, 2 * j + sj
        # P_ab = sum_n w_n <n|c_a|m><m|c^+_b|n> / (iw - E_m + E_n)
        P[i, si, j, sj, 0] = np.einsum('n,nm,mn,mnw->w', w, c[a], cdag[b],
                                       1 / (iw[None, None, :] - E[:, None, None] + E[None, :, None]))
        # H_ab = sum_n w_n <n|c^+_b|m><m|c_a|n> / (iw + E_m - E_n)
        H[i, si, j, sj, 0] = np.einsum('n,nm,mn,mnw->w', w, cdag[b], c[a],
                                       1 / (iw[None, None, :] + E[:, None, None] - E[None, :, None]))
    G = calc_one_body_green(_spectra(P, H))
    assert np.allclose(G, P + H, atol=1e-13)
    # independent check: G_ab(iw) = [(iw - H_1body)^{-1}]-like Lehmann sum over both orderings
    for (i, si), (j, sj) in itertools.product(itertools.product(range(2), range(2)), repeat=2):
        a, b = 2 * i + si, 2 * j + sj
        ref = np.einsum('n,nm,mn,mnw->w', w, c[a], cdag[b], 1 / (iw[None, None, :] - E[:, None, None] + E[None, :, None])) \
            + np.einsum('n,nm,mn,mnw->w', w, cdag[b], c[a], 1 / (iw[None, None, :] + E[:, None, None] - E[None, :, None]))
        assert np.allclose(G[i, si, j, sj, 0], ref, atol=1e-13)


def test_omp_env(monkeypatch):
    from dcore.impurity_solvers.hphi_spectrum import omp_env, available_cores
    monkeypatch.delenv('OMP_NUM_THREADS', raising=False)
    monkeypatch.delenv('OMP_WAIT_POLICY', raising=False)
    # default: available cores divided by the number of parallel HPhi processes
    env = omp_env(2)
    assert env['OMP_NUM_THREADS'] == str(max(1, available_cores() // 2))
    assert env['OMP_WAIT_POLICY'] == 'PASSIVE'
    assert omp_env(10 ** 6)['OMP_NUM_THREADS'] == '1'
    # the environment is respected
    monkeypatch.setenv('OMP_NUM_THREADS', '3')
    monkeypatch.setenv('OMP_WAIT_POLICY', 'ACTIVE')
    env = omp_env(2)
    assert env['OMP_NUM_THREADS'] == '3' and env['OMP_WAIT_POLICY'] == 'ACTIVE'
    # the solver parameter has the highest priority
    assert omp_env(2, omp_num_threads=5)['OMP_NUM_THREADS'] == '5'
