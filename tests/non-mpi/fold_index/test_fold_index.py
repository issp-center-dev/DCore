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
Unit tests for the fancy-index up/down-fold of SumkDFT_opt.

The contiguous-projection fast path (a plain slice) must produce exactly the
same result as the advanced-index (numpy.ix_) formulation, for contiguous,
non-contiguous and empty projection indices and for several `fac` values.
"""
import types

import numpy as np
import pytest

# SumkDFT_opt pulls in the backend (dcorelib/TRIQS) via dcore._dispatcher.
SumkDFT_opt = pytest.importorskip("dcore.sumkdft_opt").SumkDFT_opt


class _FakeGf:
    """Minimal stand-in for the backend Gf: only .data/.copy()/.zero()."""
    def __init__(self, data):
        self.data = data

    def copy(self):
        return _FakeGf(self.data.copy())

    def zero(self):
        self.data[...] = 0.0


def _fake_self(projindex, n_orb):
    s = types.SimpleNamespace()
    s.SO = 0
    s.spin_names_to_ind = {0: {"up": 0}}
    s.n_orbitals = np.array([[n_orb]])              # [ik, isp]
    s.corr_shells = [{"dim": len(projindex)}]
    pi = np.zeros((1, 1, 1, len(projindex)), dtype=int)
    pi[0, 0, 0, :] = projindex
    s.proj_index = pi
    return s


def _rand(shape, rng):
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


@pytest.mark.parametrize("projindex", [[1, 2, 3], [1, 3, 4], []])
@pytest.mark.parametrize("fac", [1.0, -1.0, 0.5])
def test_upfold_index_matches_advanced_index(projindex, fac):
    rng = np.random.RandomState(0)
    niw, n_orb = 4, 6
    projindex = np.array(projindex, dtype=int)
    dim = projindex.size

    inp0 = _rand((niw, n_orb, n_orb), rng)
    src = _rand((niw, dim, dim), rng)

    # reference: the original advanced-index update
    ref = inp0.copy()
    if dim > 0:
        ref[np.ix_(range(niw), projindex, projindex)] += src * fac

    got = _FakeGf(inp0.copy())
    SumkDFT_opt.upfold_index(
        _fake_self(projindex, n_orb), ik=0, ish=0, bname="up",
        gf_to_upfold=_FakeGf(src), gf_inp=got, overwrite_gf_inp=True, fac=fac)

    assert np.allclose(got.data, ref, atol=1e-13)


@pytest.mark.parametrize("projindex", [[1, 2, 3], [1, 3, 4], []])
@pytest.mark.parametrize("fac", [1.0, -1.0, 0.5])
def test_downfold_index_matches_advanced_index(projindex, fac):
    rng = np.random.RandomState(1)
    niw, n_orb = 4, 6
    projindex = np.array(projindex, dtype=int)
    dim = projindex.size

    src = _rand((niw, n_orb, n_orb), rng)
    inp0 = _rand((niw, dim, dim), rng)

    # reference: the original advanced-index update into a zeroed dim x dim block
    ref = np.zeros((niw, dim, dim), dtype=np.complex128)
    if dim > 0:
        ref[:, :, :] += src[:, projindex, :][:, :, projindex] * fac

    # overwrite_gf_inp=False returns a new (zeroed-then-filled) gf
    out = SumkDFT_opt.downfold_index(
        _fake_self(projindex, n_orb), ik=0, ish=0, bname="up",
        gf_to_downfold=_FakeGf(src), gf_inp=_FakeGf(inp0.copy()), overwrite_gf_inp=False, fac=fac)

    assert np.allclose(out.data, ref, atol=1e-13)
