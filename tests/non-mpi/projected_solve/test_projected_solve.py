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
Unit test for SumkDFT_opt._accumulate_Gloc_via_solve.

extract_G_loc only needs the correlated blocks of the lattice GF, so it solves
the lattice problem for those columns instead of forming the full inverse. This
test checks that, for a band space larger than the correlated subspace
(n_band > dim, the ab-initio downfolding regime that the existing model tests do
not cover), the solve path reproduces the full-inverse-then-downfold result.
"""
import types

import numpy as np
import pytest

SumkDFT_opt = pytest.importorskip("dcore.sumkdft_opt").SumkDFT_opt


class _FakeGf:
    def __init__(self, data):
        self.data = data


class _FakeBlockGf:
    def __init__(self, blocks):
        self._b = blocks            # dict: bname -> _FakeGf

    def __iter__(self):
        return iter(self._b.items())

    def __getitem__(self, bname):
        return self._b[bname]


def test_accumulate_Gloc_via_solve_matches_full_inverse():
    rng = np.random.RandomState(0)
    n_iw, n_band = 5, 6
    bname = "up"
    # two correlated shells living in a 6-band space; dims 2 and 1, with
    # non-contiguous band indices to also exercise that path
    shells = [{"dim": 2, "proj": [1, 3]}, {"dim": 1, "proj": [4]}]
    bz_w = 0.37
    ik = 0

    def M_data():
        A = rng.standard_normal((n_iw, n_band, n_band)) \
            + 1j * rng.standard_normal((n_iw, n_band, n_band))
        return A + n_band * np.eye(n_band)        # well conditioned

    Md = M_data()
    M = _FakeBlockGf({bname: _FakeGf(Md)})

    G_loc = [_FakeBlockGf({bname: _FakeGf(np.zeros((n_iw, s["dim"], s["dim"]), dtype=complex))})
             for s in shells]

    # fake self carrying just the attributes the method reads
    s = types.SimpleNamespace()
    s.SO = 0
    s.spin_names_to_ind = {0: {bname: 0}}
    s.n_corr_shells = len(shells)
    s.corr_shells = [{"dim": sh["dim"]} for sh in shells]
    s.bz_weights = np.array([bz_w])
    pidx_max = max(sh["dim"] for sh in shells)
    proj = np.zeros((1, 1, len(shells), pidx_max), dtype=int)
    for i, sh in enumerate(shells):
        proj[0, 0, i, 0:sh["dim"]] = sh["proj"]
    s.proj_index = proj

    SumkDFT_opt._accumulate_Gloc_via_solve(s, ik, M, G_loc)

    # reference: full inverse, then take the correlated [proj, proj] block * bz_weight
    Ginv = np.linalg.inv(Md)
    for i, sh in enumerate(shells):
        p = sh["proj"]
        ref = bz_w * Ginv[:, p, :][:, :, p]
        got = G_loc[i][bname].data
        assert np.allclose(got, ref, atol=1e-12), "shell %d mismatch" % i
