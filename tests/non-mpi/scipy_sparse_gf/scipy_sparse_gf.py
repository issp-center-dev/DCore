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
import scipy.sparse as sp
import pytest
from dcore.impurity_solvers.scipy_sparse_main import calc_gf_Lehmann


def _reference_lehmann(iws, Cdag, spin_conserve, eigvec, E_n, evx, vvx, pm):
    nf = Cdag.size; niw = iws.size; dim = evx.size
    gf = np.zeros((nf, nf, niw), dtype=complex)
    cdag_im = np.empty((nf, dim), dtype=complex)
    for i in range(nf):
        cdag_im[i] = vvx.conj().T @ Cdag[i] @ eigvec
    for l, iw in enumerate(iws):
        den = 1/(iw - evx + E_n) if pm == +1 else 1/(iw + evx - E_n)
        for i, j in np.ndindex(nf, nf):
            if spin_conserve and (i // (nf//2) != j // (nf//2)):
                continue
            gf[i, j, l] = np.einsum("m,m,m", cdag_im[i].conj(), den, cdag_im[j])
    return gf


def _rand_case(seed, nf=4, dim=6, niw=5):
    rng = np.random.default_rng(seed)
    def cx(shape):
        return rng.standard_normal(shape) + 1j*rng.standard_normal(shape)
    Cdag = np.empty(nf, dtype=object)
    for i in range(nf):
        Cdag[i] = sp.csr_matrix(cx((dim, dim)))
    eigvec = cx(dim)
    evx = rng.standard_normal(dim)
    vvx = cx((dim, dim))
    iws = 1j * (2*np.arange(niw)+1) * np.pi / 10.0
    return iws, Cdag, eigvec, evx, vvx


@pytest.mark.parametrize("pm", [+1, -1])
@pytest.mark.parametrize("spin_conserve", [True, False])
def test_lehmann_matches_reference(pm, spin_conserve):
    iws, Cdag, eigvec, evx, vvx = _rand_case(
        seed=(0 if pm == +1 else 1) + (10 if spin_conserve else 0)
    )
    ref = _reference_lehmann(iws, Cdag, spin_conserve, eigvec, 0.3, evx, vvx, pm)
    got = calc_gf_Lehmann(iws, Cdag, spin_conserve, eigvec, 0.3, evx, vvx, pm)
    assert got.shape == ref.shape
    assert np.allclose(got, ref, atol=1e-12)
