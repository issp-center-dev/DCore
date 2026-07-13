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
import os
import importlib.util
import itertools

import numpy as np

# Load hphi_spectrum.py directly by file path. The module only depends on numpy,
# so the test does not require the optional impurity-solver backends
# (dcorelib/TRIQS) that the package __init__ would otherwise pull in.
_HERE = os.path.dirname(os.path.abspath(__file__))
_MODULE_PATH = os.path.normpath(
    os.path.join(_HERE, "..", "..", "..", "src", "dcore", "impurity_solvers", "hphi_spectrum.py"))
_spec = importlib.util.spec_from_file_location("hphi_spectrum", _MODULE_PATH)
hphi_spectrum = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hphi_spectrum)
calc_one_body_green = hphi_spectrum.calc_one_body_green


def _composite(site, sigma, n_sigma=2):
    return site * n_sigma + sigma


def _build_core(G, n_site, n_sigma=2):
    """Synthesize the HPhi excitation 'core' tensor (upper triangle only) for a
    target Green's-function matrix G, following the reconstruction convention:

        B = c_i + c_j    ->  <<B|B+>> = G_ii + G_jj + G_ij + G_ji    (idx=1)
        A = c_i + i c_j  ->  <<A|A+>> = G_ii + G_jj + i(G_ji - G_ij)  (idx=0)

    The lower triangle and the diagonal idx=1 channel are left as zeros, exactly
    as the reduced parallel driver (`calc_one_body_green_core_parallel`) produces
    them.
    """
    n_flg, n_excitation = 2, 2
    n_T, n_omega = G.shape[2], G.shape[3]
    core = np.zeros(
        (n_site, n_sigma, n_site, n_sigma, n_flg, n_excitation, n_T, n_omega),
        dtype=np.complex128)
    for si, sgi, sj, sgj in itertools.product(
            range(n_site), range(n_sigma), range(n_site), range(n_sigma)):
        a, b = _composite(si, sgi, n_sigma), _composite(sj, sgj, n_sigma)
        if a > b:
            continue  # lower triangle is not computed by the reduced driver
        if a == b:
            core[si, sgi, sj, sgj, 0, 0] = G[a, a]  # diagonal uses idx=0 only
            continue
        diag = G[a, a] + G[b, b]
        core[si, sgi, sj, sgj, 1, 0] = diag + G[a, b] + G[b, a]
        core[si, sgi, sj, sgj, 0, 0] = diag + 1j * (G[b, a] - G[a, b])
    return core


def test_offdiagonal_reconstruction_fills_both_triangles():
    """G_ij AND the transposed G_ji must be recovered, even though the lower
    triangle is never computed (regression test for the dead-store bug)."""
    n_site, n_sigma, n_T, n_omega = 3, 2, 1, 4
    N = n_site * n_sigma
    rng = np.random.RandomState(2)
    G = rng.rand(N, N, n_T, n_omega) + 1j * rng.rand(N, N, n_T, n_omega)

    out = calc_one_body_green(_build_core(G, n_site, n_sigma))

    for si, sgi, sj, sgj in itertools.product(
            range(n_site), range(n_sigma), range(n_site), range(n_sigma)):
        a, b = _composite(si, sgi, n_sigma), _composite(sj, sgj, n_sigma)
        if a == b:
            continue
        assert np.allclose(out[si, sgi, sj, sgj], G[a, b], atol=1e-12), \
            "off-diagonal mismatch at (%d,%d,%d,%d)" % (si, sgi, sj, sgj)


def test_diagonal_reconstruction():
    n_site, n_sigma, n_T, n_omega = 2, 2, 1, 3
    N = n_site * n_sigma
    rng = np.random.RandomState(5)
    G = rng.rand(N, N, n_T, n_omega) + 1j * rng.rand(N, N, n_T, n_omega)

    out = calc_one_body_green(_build_core(G, n_site, n_sigma))

    for si, sgi in itertools.product(range(n_site), range(n_sigma)):
        a = _composite(si, sgi, n_sigma)
        assert np.allclose(out[si, sgi, si, sgi], G[a, a], atol=1e-12)


class _SerialExecutor:
    """Drop-in stand-in for ProcessPoolExecutor that runs map() in-process, so
    the monkeypatched worker is actually used (no subprocess / pickling)."""
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def map(self, fn, iterable):
        return [fn(x) for x in iterable]


def test_parallel_driver_prunes_tasks_and_reconstructs(monkeypatch):
    """Drive calc_one_body_green_core_parallel end-to-end with a mocked HPhi
    worker: assert only the upper triangle is scheduled (diagonal idx=1 skipped)
    and the assembled Green's function matches the target."""
    import concurrent.futures
    import shutil

    n_site, n_sigma, n_T, n_omega = 3, 2, 1, 4
    N = n_site * n_sigma
    rng = np.random.RandomState(7)
    G = rng.rand(N, N, n_T, n_omega) + 1j * rng.rand(N, N, n_T, n_omega)

    scheduled = []

    def fake_worker(task):
        sitei, sigmai, sitej, sigmaj, idx, ex_state, _p = task
        scheduled.append((sitei, sigmai, sitej, sigmaj, idx, ex_state))
        a = _composite(sitei, sigmai, n_sigma)
        b = _composite(sitej, sigmaj, n_sigma)
        out = np.zeros((n_T, n_omega), dtype=np.complex128)
        if a == b:
            if idx == 0 and ex_state == 0:
                out[:] = G[a, a]
            return out
        if ex_state == 0:
            diag = G[a, a] + G[b, b]
            if idx == 1:
                out[:] = diag + G[a, b] + G[b, a]
            else:
                out[:] = diag + 1j * (G[b, a] - G[a, b])
        return out

    monkeypatch.setattr(hphi_spectrum, "check_eta", lambda p_common: None)
    monkeypatch.setattr(hphi_spectrum, "calc_one_body_green_core", fake_worker)
    monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", _SerialExecutor)
    monkeypatch.setattr(shutil, "rmtree", lambda *a, **k: None)

    p_common = (n_site, [1.0], 1, 1e-4, "./HPhi", "zvo", "./output", 1)
    out = hphi_spectrum.calc_one_body_green_core_parallel(p_common, max_workers=1)

    # Pruning: no lower-triangle pairs, and no diagonal idx==1 tasks.
    for sitei, sigmai, sitej, sigmaj, idx, ex_state in scheduled:
        a = _composite(sitei, sigmai, n_sigma)
        b = _composite(sitej, sigmaj, n_sigma)
        assert a <= b, "lower-triangle task scheduled"
        if a == b:
            assert idx == 0, "diagonal idx==1 task should be skipped"

    # Reconstruction over the full off-diagonal block (including the never-computed
    # lower triangle).
    for si, sgi, sj, sgj in itertools.product(
            range(n_site), range(n_sigma), range(n_site), range(n_sigma)):
        a, b = _composite(si, sgi, n_sigma), _composite(sj, sgj, n_sigma)
        if a == b:
            continue
        assert np.allclose(out[si, sgi, sj, sgj], G[a, b], atol=1e-12)


def test_parallel_driver_rejects_empty_tasks(monkeypatch):
    monkeypatch.setattr(hphi_spectrum, "check_eta", lambda p_common: None)
    p_common = (0, [1.0], 1, 1e-4, "./HPhi", "zvo", "./output", 1)
    try:
        hphi_spectrum.calc_one_body_green_core_parallel(p_common, max_workers=1)
    except ValueError:
        return
    raise AssertionError("expected ValueError for n_site == 0")
