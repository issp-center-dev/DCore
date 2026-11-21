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
Pure-Python checks for the Stage-3 direct bra/ket launch enumeration
(``_braket_sector_jobs``), CI-safe (no HPhi executable needed).

Each launch is one (ex_state, spin) sector with kets = bras = all same-spin
sites; cross-spin operators are never emitted, and a canonical occupancy drops
the (ex_state, spin) sectors whose excited state vanishes.
"""

import os

import numpy as np
import pytest

from dcore.impurity_solvers.hphi_spectrum import (
    _braket_sector_jobs, _braket_ne_sector_jobs, _braket_store_op, _braket_assemble,
    _braket_one_body_green, CalcSpectrumCore, calc_one_body_green_core_parallel)


def test_grand_canonical_runs_every_exstate_and_spin():
    n_site = 3
    jobs = _braket_sector_jobs(n_site, None)
    # 2 ex_state x 2 spin = 4 launches.
    assert len(jobs) == 4
    seen = set()
    for ex_state, ket_ops, bra_ops in jobs:
        assert ket_ops == bra_ops, "ket and bra operator lists must match"
        spins = {sg for (_, sg) in ket_ops}
        assert len(spins) == 1, "a launch is a single spin sector"
        sigma = spins.pop()
        # all sites of that spin are present
        assert [s for (s, _) in ket_ops] == list(range(n_site))
        seen.add((ex_state, sigma))
    assert seen == {(e, s) for e in (0, 1) for s in (0, 1)}


def test_no_cross_spin_operators():
    for occ in (None, (1, 1, 2), (2, 0, 2)):
        for ex_state, ket_ops, bra_ops in _braket_sector_jobs(2, occ):
            sigmas = {sg for (_, sg) in ket_ops} | {sg for (_, sg) in bra_ops}
            assert len(sigmas) == 1, "cross-spin operators must never be batched"


def test_canonical_drops_empty_and_full_spin_sectors():
    n_site = 2
    # sector with up fully empty (n_up=0) and down half-filled (n_down=1)
    jobs = _braket_sector_jobs(n_site, (0, 1, n_site))
    got = {(ex_state, ket_ops[0][1]) for ex_state, ket_ops, _ in jobs}
    # up (sigma 0): annihilation (ex 0) invalid (empty); creation (ex 1) valid
    assert (0, 0) not in got
    assert (1, 0) in got
    # down (sigma 1, n_down=1): both annihilation and creation valid (0 < 1 < 2)
    assert (0, 1) in got and (1, 1) in got

    # full up sector (n_up = n_site): creation invalid, annihilation valid
    jobs_full = _braket_sector_jobs(n_site, (n_site, 1, n_site))
    got_full = {(ex_state, ket_ops[0][1]) for ex_state, ket_ops, _ in jobs_full}
    assert (1, 0) not in got_full   # up creation on a full spin vanishes
    assert (0, 0) in got_full       # up annihilation valid


def test_braket_reader_handles_single_site_unsuffixed_files(tmp_path):
    """n_site == 1 gives num_op == num_bra == 1, so HPhi writes the unsuffixed
    <header>_DynamicalGreen_<idx>.dat; the braket reader must accept that layout."""
    core = CalcSpectrumCore([1.0], 2, 1e-4, header="zvo", output_dir="output")
    core.energy_list = [0.0, 0.5]
    core.ene_min = 0.0
    spec_dir = os.path.join(str(tmp_path), "output")
    os.makedirs(spec_dir)
    # 3 frequencies; columns: re(w) im(w) re(G) im(G). idx encodes the value so we can check it.
    for idx in range(2):
        rows = np.array([[w, 0.0, float(idx) + w, -float(idx)] for w in range(3)])
        np.savetxt(os.path.join(spec_dir, "zvo_DynamicalGreen_{}.dat".format(idx)), rows)
    freqs, finite = core._finite_T_spectrum_braket(str(tmp_path), num_op=1, num_bra=1)
    assert (0, 0, 1.0) in finite          # the (op=0, bra=0) diagonal element exists
    assert finite[(0, 0, 1.0)].shape == (3,)
    assert len(freqs) == 3


def test_ne_sector_jobs_use_both_spins_and_include_cross_spin():
    """The Ne-only (spin-orbit) route batches BOTH spins per launch, so cross-spin operator
    pairs are present (unlike the same-spin (Ne, 2Sz) route)."""
    n_site = 2
    jobs = _braket_ne_sector_jobs(n_site, ne=2)  # interior Ne: both ex_state valid
    assert len(jobs) == 2
    for ex_state, ket_ops, bra_ops in jobs:
        assert ket_ops == bra_ops
        # all 2*n_site spin-orbitals present, both spins
        assert sorted(ket_ops) == [(s, sg) for s in range(n_site) for sg in (0, 1)]
        assert {sg for (_, sg) in ket_ops} == {0, 1}


def test_ne_sector_jobs_drop_empty_and_full():
    n_site = 2
    n_so = 2 * n_site
    # vacuum (Ne=0): only creation (ex_state 1)
    assert [e for e, _, _ in _braket_ne_sector_jobs(n_site, 0)] == [1]
    # full (Ne = 2*n_site): only annihilation (ex_state 0)
    assert [e for e, _, _ in _braket_ne_sector_jobs(n_site, n_so)] == [0]
    # interior: both
    assert sorted(e for e, _, _ in _braket_ne_sector_jobs(n_site, 2)) == [0, 1]


def test_ne_sector_jobs_reject_out_of_range_ne():
    with pytest.raises(ValueError, match=r"ne must be in"):
        _braket_ne_sector_jobs(2, ne=-1)
    with pytest.raises(ValueError, match=r"ne must be in"):
        _braket_ne_sector_jobs(2, ne=5)  # 2*n_site = 4
    # n_site = 0: only the empty (Ne=0) sector; creation valid (Ne <= 2*n_site-1 = -1 is False),
    # annihilation invalid -> no jobs.
    assert _braket_ne_sector_jobs(0, ne=0) == []


def test_braket_store_op_convention():
    # g = G_{bra, ket}. Annihilation (ex 0): stored at [ket][bra]. Creation (ex 1): transposed [bra][ket].
    bra, ket = (0, 0), (1, 0)  # same spin, off-diagonal
    assert _braket_store_op(bra, ket, 0) == (ket, bra)   # (1,0),(0,0)
    assert _braket_store_op(bra, ket, 1) == (bra, ket)   # (0,0),(1,0)  -> transpose
    # diagonal: both channels store at the same place
    d = (0, 0)
    assert _braket_store_op(d, d, 0) == _braket_store_op(d, d, 1) == (d, d)
    # cross-spin uses the identical convention
    bra_x, ket_x = (0, 0), (0, 1)
    assert _braket_store_op(bra_x, ket_x, 0) == (ket_x, bra_x)
    assert _braket_store_op(bra_x, ket_x, 1) == (bra_x, ket_x)


def test_braket_assemble_places_same_and_cross_spin(monkeypatch):
    """Mocked accumulation: feed _braket_assemble per-launch {(row_op,col_op): g} dicts (the
    storage keys _braket_store_op produces) for both an annihilation and a creation launch over
    BOTH spins, and check every (i,si,j,sj) slot -- including cross-spin -- lands where expected.
    g carries a unique tag per (bra, ket, ex_state) so mis-routing would be caught."""
    n_site = 2
    ops = [(s, sg) for s in range(n_site) for sg in range(2)]  # all 4 spin-orbitals
    # tag(bra, ket, ex) -> a distinct (n_T=1, n_omega=1) complex "spectrum"
    def tag(bra, ket, ex):
        code = (bra[0] * 2 + bra[1]) * 10 + (ket[0] * 2 + ket[1]) + 100 * ex
        return np.array([[code + 0j]])

    results = []
    expected = {}  # (row_op, col_op) -> summed g
    for ex in (0, 1):
        res = {}
        for ket in ops:
            for bra in ops:
                g = tag(bra, ket, ex)
                key = _braket_store_op(bra, ket, ex)  # storage index
                res[key] = g
                expected[key] = expected.get(key, np.zeros((1, 1), dtype=complex)) + g
        results.append(res)

    obg = _braket_assemble(results, n_site)
    assert obg.shape == (n_site, 2, n_site, 2, 1, 1)
    # every stored slot matches the summed expectation, incl. cross-spin (si != sj)
    saw_cross = False
    for (row_op, col_op), val in expected.items():
        (ri, rsg), (ci, csg) = row_op, col_op
        assert obg[ri][rsg][ci][csg][0, 0] == val[0, 0]
        if rsg != csg:
            saw_cross = True
    assert saw_cross, "cross-spin slots must be exercised"


def test_braket_one_body_green_rejects_ne_only_with_occupancy():
    # ne_only and sector_occupancy are mutually exclusive; the guard runs before any HPhi launch.
    with pytest.raises(ValueError, match="mutually exclusive"):
        _braket_one_body_green(_p_common(2), max_workers=1,
                               sector_occupancy=(1, 1, 2), ne_only=2)


def _p_common(n_site):
    return (n_site, [1.0], 1, 1e-4, "./HPhi", "zvo", "./output", 1)


def test_braket_rejects_grand_canonical(monkeypatch):
    """DCORE_HPHI_BRAKET=1 with no sector (grand canonical) must fail cleanly, not run."""
    monkeypatch.setenv("DCORE_HPHI_BRAKET", "1")
    with pytest.raises(RuntimeError, match="grand-canonical bra/ket route is unsupported"):
        calc_one_body_green_core_parallel(_p_common(2), max_workers=1, sector_occupancy=None)


def test_braket_rejects_zero_site(monkeypatch):
    """n_site == 0 raises the stable ValueError on the bra/ket path too (not StopIteration)."""
    monkeypatch.setenv("DCORE_HPHI_BRAKET", "1")
    with pytest.raises(ValueError, match="n_site must be >= 1"):
        calc_one_body_green_core_parallel(_p_common(0), max_workers=1, sector_occupancy=(0, 0, 0))
