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
import pytest
sparse_ir = pytest.importorskip("sparse_ir")

from dcore.ir_basis import get_basis


def test_get_basis_builds():
    b = get_basis(beta=10.0, wmax=20.0, eps=1e-10, statistics='F')
    assert isinstance(b, sparse_ir.FiniteTempBasis)
    assert b.size > 0


def test_get_basis_is_cached():
    b1 = get_basis(beta=10.0, wmax=20.0, eps=1e-10, statistics='F')
    b2 = get_basis(beta=10.0, wmax=20.0, eps=1e-10, statistics='F')
    assert b1 is b2  # same object returned from cache


def test_get_basis_distinct_keys():
    b1 = get_basis(beta=10.0, wmax=20.0, eps=1e-10, statistics='F')
    b3 = get_basis(beta=10.0, wmax=30.0, eps=1e-10, statistics='F')
    assert b1 is not b3
