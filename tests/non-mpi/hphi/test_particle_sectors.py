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
"""Unit tests for the (Ne, 2Sz) particle-number sector enumeration used by the
canonical (sector-resolved) HPhi spectrum path."""

from math import comb

import pytest

from dcore.impurity_solvers.hphi import enumerate_particle_sectors


@pytest.mark.parametrize("n_site", [1, 2, 3, 4])
def test_sectors_partition_the_grand_canonical_space(n_site):
    sectors = enumerate_particle_sectors(n_site)
    # the sector dimensions must add up to the full grand-canonical space 4**n_site
    assert sum(s['dim'] for s in sectors) == 4 ** n_site
    # (Ne, 2Sz) labels are unique, and there are exactly (n_site+1)**2 of them
    labels = [(s['Ne'], s['two_Sz']) for s in sectors]
    assert len(labels) == len(set(labels)) == (n_site + 1) ** 2
    # the list is returned sorted by (Ne, 2Sz) -- the per-sector loop relies on this order
    assert labels == sorted(labels)


def test_n_site_zero_is_the_single_vacuum_sector():
    assert enumerate_particle_sectors(0) == [
        {'Ne': 0, 'two_Sz': 0, 'n_up': 0, 'n_down': 0, 'dim': 1}
    ]


def test_negative_n_site_is_rejected():
    with pytest.raises(ValueError):
        enumerate_particle_sectors(-1)


def test_sector_quantum_numbers_are_consistent():
    for s in enumerate_particle_sectors(3):
        assert s['n_up'] + s['n_down'] == s['Ne']
        assert s['n_up'] - s['n_down'] == s['two_Sz']
        assert 0 <= s['n_up'] <= 3 and 0 <= s['n_down'] <= 3
        assert s['dim'] == comb(3, s['n_up']) * comb(3, s['n_down'])


def test_half_filled_sz0_is_the_largest_sector():
    # for EVEN n_site the (Ne = n_site, 2Sz = 0) sector is the largest (half filling, Sz=0);
    # for odd n_site, Ne and 2Sz must share parity so 2Sz=0 is not allowed at Ne=n_site.
    n_site = 4
    sectors = enumerate_particle_sectors(n_site)
    largest = max(sectors, key=lambda s: s['dim'])
    assert largest['Ne'] == n_site
    assert largest['two_Sz'] == 0
