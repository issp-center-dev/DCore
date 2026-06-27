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
Unit tests for the HPhi off-diagonal excitation definition (no HPhi binary).

These guard the fix to the complex excitation A = c_i + i c_j: the creation
channel (ex_state == 1) must use the conjugate operator A^dag = c_i^dag - i c_j^dag,
i.e. the imaginary coefficient of c_j flips sign between the annihilation and
creation channels. Getting this wrong gives the anti-symmetric (imaginary) part
of the off-diagonal Green's function a divergent high-frequency tail.
"""

import os

import pytest

from dcore.impurity_solvers.hphi_spectrum import CalcSpectrumCore


def _parse_single_ex(path):
    """Return the operator rows (site, spin, op, re, im) from a single_ex.def file."""
    rows = []
    with open(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) == 5:  # site spin op re im
                rows.append((int(tok[0]), int(tok[1]), int(tok[2]),
                             float(tok[3]), float(tok[4])))
    return rows


def _write_offdiag_excitation(tmp_path, ex_state):
    core = CalcSpectrumCore([0.1], 1, 1e-4, path_to_HPhi="HPhi")
    # off-diagonal pair: (site 0, spin 0) and (site 1, spin 0) -> composites 0 != 2
    core._make_single_excitation(0, 0, 1, 0, file_name="single_ex.def",
                                 ex_state=ex_state, flg_complex=True,
                                 calc_dir=str(tmp_path))
    return _parse_single_ex(os.path.join(str(tmp_path), "single_ex.def"))


def test_annihilation_channel_uses_plus_i(tmp_path):
    rows = _write_offdiag_excitation(tmp_path, ex_state=0)
    assert len(rows) == 2
    # c_i with coefficient 1, c_j with coefficient +i
    assert rows[0][3:] == (1.0, 0.0)
    assert rows[1][3:] == (0.0, 1.0)
    assert all(op == 0 for _, _, op, _, _ in rows)


def test_creation_channel_uses_minus_i(tmp_path):
    rows = _write_offdiag_excitation(tmp_path, ex_state=1)
    assert len(rows) == 2
    # A^dag: c_i^dag with coefficient 1, c_j^dag with coefficient -i
    assert rows[0][3:] == (1.0, 0.0)
    assert rows[1][3:] == (0.0, -1.0)
    assert all(op == 1 for _, _, op, _, _ in rows)


@pytest.mark.parametrize("ex_state", [0, 1])
def test_real_combination_is_unaffected(tmp_path, ex_state):
    # flg_complex=False -> B = c_i + c_j, both coefficients real and independent of channel
    core = CalcSpectrumCore([0.1], 1, 1e-4, path_to_HPhi="HPhi")
    core._make_single_excitation(0, 0, 1, 0, file_name="b.def",
                                 ex_state=ex_state, flg_complex=False,
                                 calc_dir=str(tmp_path))
    rows = _parse_single_ex(os.path.join(str(tmp_path), "b.def"))
    assert [r[3:] for r in rows] == [(1.0, 0.0), (1.0, 0.0)]
