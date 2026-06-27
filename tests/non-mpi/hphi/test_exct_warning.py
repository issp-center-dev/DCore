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
Unit tests for the HPhi 'exct too small' warning (pure Python, no HPhi binary).
"""

import numpy
import pytest

from dcore.impurity_solvers.hphi import (
    read_eigenenergies,
    warn_if_exct_truncates_thermal_trace,
)

# Exact 2-orbital Kanamori atomic spectrum (eps=-0.5, U=4, U'=0.8, J=0): 16 states.
KANAMORI_SPECTRUM = numpy.array(
    [-0.5] * 4 + [-0.2] * 4 + [0.0] + [3.0] * 2 + [4.1] * 4 + [9.2]
)
EXCT_MAX = 16
BETA = 10.0


def _write_energy_file(path, energies):
    with open(path, "w") as f:
        for i, e in enumerate(numpy.sort(energies)):
            f.write(f"State {i}\n  Energy  {e:.16f} \n  Doublon  0.0 \n  Sz  0.0 \n")


def test_read_eigenenergies(tmp_path):
    fn = str(tmp_path / "zvo_energy.dat")
    _write_energy_file(fn, KANAMORI_SPECTRUM[:5])
    e = read_eigenenergies(fn)
    assert e.size == 5
    assert numpy.isclose(e.min(), -0.5)


def _warns(tmp_path, exct):
    fn = str(tmp_path / "zvo_energy.dat")
    _write_energy_file(fn, KANAMORI_SPECTRUM[: min(exct, EXCT_MAX)])
    return fn


@pytest.mark.parametrize("exct", [1, 4, 9])
def test_warns_when_exct_too_small(tmp_path, capsys, exct):
    fn = _warns(tmp_path, exct)
    warn_if_exct_truncates_thermal_trace(fn, BETA, exct, EXCT_MAX)
    err = capsys.readouterr().err
    assert "exct" in err and "may be too small" in err


@pytest.mark.parametrize("exct", [11, 16])
def test_silent_when_exct_sufficient(tmp_path, capsys, exct):
    fn = _warns(tmp_path, exct)
    warn_if_exct_truncates_thermal_trace(fn, BETA, exct, EXCT_MAX)
    err = capsys.readouterr().err
    assert err == ""


def test_no_warning_when_full_space(tmp_path, capsys):
    # exct == exct_max must never warn, even if the last state has weight.
    fn = _warns(tmp_path, EXCT_MAX)
    warn_if_exct_truncates_thermal_trace(fn, BETA, EXCT_MAX, EXCT_MAX)
    assert capsys.readouterr().err == ""


def test_missing_file_does_not_raise(tmp_path, capsys):
    # A missing energy file must not abort the solver.
    warn_if_exct_truncates_thermal_trace(
        str(tmp_path / "does_not_exist.dat"), BETA, 1, EXCT_MAX
    )
    assert capsys.readouterr().err == ""
