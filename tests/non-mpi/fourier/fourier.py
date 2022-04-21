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

from dcore.fourier import _fft_fermion_w2t, _fft_fermion_t2w, _matsubara_freq_fermion, bgf_fourier_w2t
from dcore.tools import make_block_gf
from dcore._dispatcher import BlockGf, Gf, GfImFreq, GfImTime
import numpy
import os
from itertools import product


def _make_g0():
    beta = 10.0
    nw = 1024
    nt = 2 * nw
    e0 = 2.0
    a = 0.1

    iw = _matsubara_freq_fermion(beta, nw)
    g0_w = a / (iw - e0)

    t = numpy.linspace(0, beta, nt+1)
    g0_t = -numpy.exp(-t*e0) / (numpy.exp(-beta*e0) + 1) * a

    return g0_w, g0_t, beta, a


# test for t2w
def test_fft_fermion_t2w(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    g0_w, g0_t, beta, a = _make_g0()
    g0_w_fft = _fft_fermion_t2w(g0_t, beta)
    assert numpy.allclose(g0_w_fft, g0_w, atol=1e-4)


# test for w2t
def test_fft_fermion_w2t(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    g0_w, g0_t, beta, a = _make_g0()
    g0_t_fft = _fft_fermion_w2t(g0_w, beta, a=a)
    assert numpy.allclose(g0_t_fft, g0_t, atol=1e-4)


# test for w2t
def test_fft_bgf_w2t(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    g0_w, g0_t, beta, a = _make_g0()

    nt = g0_t.size
    nw = g0_w.size // 2

    # Set BlockGf for G(iw)
    gf_struct = {'up': [0, 1]}
    tail = {'up': numpy.identity(2) * a}

    bgf_w = make_block_gf(GfImFreq, gf_struct, beta, nw)
    for name, gf in bgf_w:
        nw_2, norb1, norb2 = gf.data.shape
        assert nw_2 == nw*2
        for i, j in product(range(norb1), range(norb2)):
            if i == j:
                gf.data[:, i, j] = g0_w[:]
            else:
                gf.data[:, i, j] = numpy.zeros(2*nw)

    # FFT
    bgf_t = bgf_fourier_w2t(bgf_w, tail=tail)

    # Compare
    for name, gf in bgf_t:
        nt_2, norb1, norb2 = gf.data.shape
        assert nt_2 == nt
        for i, j in product(range(norb1), range(norb2)):
            if i == j:
                assert numpy.allclose(gf.data[:, i, j], g0_t, atol=1e-4)
            else:
                assert numpy.allclose(gf.data[:, i, j], numpy.zeros(nt))
