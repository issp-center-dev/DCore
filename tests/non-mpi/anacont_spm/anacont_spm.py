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

import numpy as np

def test_find_sum_rule_const():
    from dcore.anacont_spm import _find_sum_rule_const
    a = 1.123
    b = 2.234
    beta = 40
    n_matsubara = 1000
    wn = np.pi / beta * (2 * np.arange(0, n_matsubara, 1) + 1)
    gf_wn = a - 1j * b / wn
    a_test, b_test = _find_sum_rule_const(wn, gf_wn, 100, False)
    assert np.allclose(a, a_test, atol=1e-9)
    assert np.allclose(b, b_test, atol=1e-9)

test_find_sum_rule_const()
