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


def test_overlapping_shells_share_one_solve(monkeypatch):
    rng = np.random.RandomState(12)
    md = rng.randn(4, 6, 6) + 1j * rng.randn(4, 6, 6) + 6 * np.eye(6)
    indices = [[3, 1], [1, 4], [4, 3]]
    s = types.SimpleNamespace(
        SO=0, spin_names_to_ind={0: {'up': 0}}, n_corr_shells=3,
        corr_shells=[{'dim': 2}] * 3, bz_weights=[0.7],
        proj_index=np.array(indices)[None, None, :, :])
    gloc = [_FakeBlockGf({'up': _FakeGf(np.ones((4, 2, 2), complex))})
            for _ in indices]
    original = np.linalg.solve
    calls = []

    def solve(a, b):
        calls.append(b.shape)
        return original(a, b)

    monkeypatch.setattr(np.linalg, 'solve', solve)
    SumkDFT_opt._accumulate_Gloc_via_solve(
        s, 0, _FakeBlockGf({'up': _FakeGf(md)}), gloc)
    inv = np.linalg.inv(md)
    for g, p in zip(gloc, indices):
        np.testing.assert_allclose(g['up'].data, 1 + 0.7 * inv[:, p, :][:, :, p])
    assert calls == [(4, 6, 3)]


def test_lattice_denominator_has_independent_frequency_cache(monkeypatch):
    # Exercise the real lattice_gf control flow with small NumPy-backed blocks.
    # A denominator at another k point must not resize the inverse's cache.
    globals_ = SumkDFT_opt.lattice_gf.__globals__

    class Gf(_FakeGf):
        def __init__(self, indices, mesh):
            super().__init__(np.zeros((2, len(indices), len(indices)), complex))
            self.mesh = mesh

    class Block(_FakeBlockGf):
        def __init__(self, name_list, block_list, make_copies):
            super().__init__(dict(zip(name_list, block_list)))
            self.mesh = block_list[0].mesh

        def zero(self):
            for _, gf in self:
                gf.data[...] = 0

        def __lshift__(self, value):
            for name, gf in self:
                gf.data[...] = (value[name].data if isinstance(value, Block)
                                else value * np.eye(gf.data.shape[1]))
            return self

        def __isub__(self, matrices):
            for (_, gf), matrix in zip(self, matrices):
                gf.data[...] -= matrix
            return self

        def invert(self):
            for _, gf in self:
                gf.data[...] = np.linalg.inv(gf.data)

    monkeypatch.setitem(globals_, 'MeshImFreq', lambda **kw: types.SimpleNamespace(beta=kw['beta']))
    monkeypatch.setitem(globals_, 'GfImFreq', Gf)
    monkeypatch.setitem(globals_, 'BlockGf', Block)
    monkeypatch.setitem(globals_, 'iOmega_n', 1j)
    s = types.SimpleNamespace(
        SO=0, spin_names_to_ind={0: {'up': 0}}, spin_block_names={0: ['up']},
        n_spin_blocks={0: 1}, n_orbitals=np.array([[2], [3]]), h_field=0,
        hopping_part=[np.zeros((1, 2, 2)), np.zeros((1, 3, 3))])
    first = SumkDFT_opt.lattice_gf(s, 0, mu=0, with_Sigma=False)._b['up'].data.copy()
    SumkDFT_opt.lattice_gf(s, 1, mu=0, with_Sigma=False, invert=False)
    again = SumkDFT_opt.lattice_gf(s, 0, mu=0, with_Sigma=False)
    np.testing.assert_allclose(again['up'].data, first)
