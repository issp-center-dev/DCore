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
    _braket_sector_jobs, CalcSpectrumCore, calc_one_body_green_core_parallel)


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


def _p_common(n_site):
    return (n_site, [1.0], 1, 1e-4, "./HPhi", "zvo", "./output", 1)


def test_braket_rejects_grand_canonical(monkeypatch):
    """DCORE_HPHI_BRAKET=1 with no sector (grand canonical) must fail cleanly, not run."""
    monkeypatch.setenv("DCORE_HPHI_BRAKET", "1")
    with pytest.raises(RuntimeError, match="canonical sector"):
        calc_one_body_green_core_parallel(_p_common(2), max_workers=1, sector_occupancy=None)


def test_braket_rejects_zero_site(monkeypatch):
    """n_site == 0 raises the stable ValueError on the bra/ket path too (not StopIteration)."""
    monkeypatch.setenv("DCORE_HPHI_BRAKET", "1")
    with pytest.raises(ValueError, match="n_site must be >= 1"):
        calc_one_body_green_core_parallel(_p_common(0), max_workers=1, sector_occupancy=(0, 0, 0))
