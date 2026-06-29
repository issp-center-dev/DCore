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
CI-safe checks for ``enumerate_ne_sectors`` -- the Ne-only (2Sz-free) sector
decomposition used by the spin-orbit-capable direct bra/ket path
(HubbardNConserved). Pure counting; no HPhi executable needed.
"""

from math import comb

import pytest

from dcore.impurity_solvers.hphi import enumerate_ne_sectors, enumerate_particle_sectors


@pytest.mark.parametrize("n_site", [1, 2, 3, 4])
def test_ne_sectors_partition_the_full_fock_space(n_site):
    sectors = enumerate_ne_sectors(n_site)
    # one sector per Ne = 0 .. 2*n_site
    assert [s['Ne'] for s in sectors] == list(range(2 * n_site + 1))
    # dimensions sum to the full grand-canonical space
    assert sum(s['dim'] for s in sectors) == 4 ** n_site
    # each Ne sector dimension is C(2*n_site, Ne)
    for s in sectors:
        assert s['dim'] == comb(2 * n_site, s['Ne'])


@pytest.mark.parametrize("n_site", [1, 2, 3, 4])
def test_ne_sector_is_union_of_2sz_sectors(n_site):
    """An Ne-only sector merges all (Ne, 2Sz) sectors of the same Ne, so its dimension
    equals the sum of the corresponding (Ne, 2Sz) sector dimensions."""
    fine = enumerate_particle_sectors(n_site)
    fine_dim_by_ne = {}
    for s in fine:
        fine_dim_by_ne[s['Ne']] = fine_dim_by_ne.get(s['Ne'], 0) + s['dim']
    for s in enumerate_ne_sectors(n_site):
        assert s['dim'] == fine_dim_by_ne[s['Ne']]


def test_n_site_zero_is_the_single_vacuum_sector():
    assert enumerate_ne_sectors(0) == [{'Ne': 0, 'dim': 1}]


def test_negative_n_site_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        enumerate_ne_sectors(-1)
