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
Unit tests for the HPhi Gf two-level parallel layout (pure Python, no HPhi).
"""

import pytest

from dcore.impurity_solvers.hphi import gf_parallel_layout

CMD = "mpirun -np 16"


def test_default_is_serial_backward_compatible():
    # n_procs_per_hphi=1 reproduces the previous behaviour: serial HPhi, np workers.
    n_inner, n_outer, prefix = gf_parallel_layout(CMD, 16, 1)
    assert (n_inner, n_outer, prefix) == (1, 16, "")


def test_unknown_np_falls_back_to_serial():
    n_inner, n_outer, prefix = gf_parallel_layout(CMD, None, 4)
    assert (n_inner, n_outer, prefix) == (1, None, "")


@pytest.mark.parametrize(
    "np_total,req,exp_inner,exp_outer",
    [
        (16, 4, 4, 4),    # clean split
        (16, 16, 16, 1),  # all ranks to one HPhi
        (64, 4, 4, 16),
    ],
)
def test_clean_power_of_four_split(np_total, req, exp_inner, exp_outer):
    n_inner, n_outer, prefix = gf_parallel_layout(CMD, np_total, req)
    assert n_inner == exp_inner
    assert n_outer == exp_outer
    assert prefix == f"mpirun -np {exp_inner}"


def test_non_power_of_four_is_rounded_down():
    # 8 is not a power of four -> rounded down to 4.
    n_inner, n_outer, prefix = gf_parallel_layout(CMD, 16, 8)
    assert n_inner == 4
    assert n_outer == 4
    assert prefix == "mpirun -np 4"


def test_inner_cannot_exceed_total():
    # request more inner ranks than available -> clamped to np_total.
    n_inner, n_outer, prefix = gf_parallel_layout(CMD, 4, 16)
    assert n_inner == 4
    assert n_outer == 1
    assert prefix == "mpirun -np 4"


def test_prefix_uses_last_token_of_command():
    # the launcher template is preserved, only the rank count is replaced.
    _, _, prefix = gf_parallel_layout("mpiexec --bind-to core -n 16", 16, 4)
    assert prefix == "mpiexec --bind-to core -n 4"
