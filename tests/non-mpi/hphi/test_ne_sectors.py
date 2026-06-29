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

from dcore.impurity_solvers.hphi import (
    enumerate_ne_sectors, enumerate_particle_sectors, ne_initial_sectors)
from dcore.impurity_solvers.hphi_spectrum import _braket_ne_sector_jobs


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


@pytest.mark.parametrize("n_site", [1, 2, 3])
def test_ne_only_route_drops_exactly_the_two_boundary_initial_sectors(n_site):
    """The spin-orbit (HubbardNConserved) bra/ket route cannot run Ne=0 (HPhi rejects Ncond=0)
    nor Ne=2*n_site (an sz() Hilbert-construction edge case) as thermally occupied INITIAL
    sectors. Lock the documented guarantee by exercising the *production* boundary filter
    ``ne_initial_sectors`` (not a local restatement): it drops exactly those two and keeps every
    interior sector."""
    n_so = 2 * n_site
    exct = 4 ** n_site  # large enough that the exct/dim gate keeps every non-empty sector
    kept, dropped = ne_initial_sectors(n_site, exct)
    assert dropped == [0, n_so]
    assert [s['Ne'] for s in kept] == list(range(1, n_so))  # every interior sector retained


@pytest.mark.parametrize("n_site", [1, 2, 3])
def test_interior_sectors_still_reach_the_boundary_excited_spaces(n_site):
    """Although Ne=0 / Ne=2*n_site are not run as initial sectors, the retained interior sectors
    still reach those boundary EXCITED spaces through their Ne+-1 transitions: Ne=1 annihilates
    into Ne=0 (ex_state 0), and Ne=2*n_site-1 creates into Ne=2*n_site (ex_state 1). This is why
    the boundary states are not simply absent from the spectrum."""
    n_so = 2 * n_site
    # Ne=1: annihilation (ex_state 0) reaches the Ne=0 boundary excited space
    ex_lo = [ex for (ex, _ket, _bra) in _braket_ne_sector_jobs(n_site, 1)]
    assert 0 in ex_lo
    # Ne=2*n_site-1: creation (ex_state 1) reaches the Ne=2*n_site boundary excited space
    ex_hi = [ex for (ex, _ket, _bra) in _braket_ne_sector_jobs(n_site, n_so - 1)]
    assert 1 in ex_hi


def test_n_site_zero_is_the_single_vacuum_sector():
    assert enumerate_ne_sectors(0) == [{'Ne': 0, 'dim': 1}]


def test_negative_n_site_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        enumerate_ne_sectors(-1)
