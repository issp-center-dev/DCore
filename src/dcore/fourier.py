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

import numpy
from scipy import fft
from itertools import product
from dcore._dispatcher import BlockGf, Gf, GfImFreq, GfImTime, MeshImFreq
from dcore.tools import make_block_gf


def _matsubara_freq_fermion(beta, nw):
    iw_positive = numpy.array([(2*i+1) * numpy.pi / beta for i in range(nw)])
    iw = numpy.append(-iw_positive[::-1], iw_positive) * 1j
    return iw


def _fft_fermion_w2t(gw, beta, a=1):
    """FFT from G(iw) to G(tau)

    Args:
        gw (numpy.ndarray(2*nw)): G(iw) including w>0 and w<0
        beta (float): Inverse temperature
        a (float, optional): Coefficient of 1/iw. Defaults to 1.

    Returns:
        numpy.ndarray(nt+1): G(tau), nt=2*nw
    """
    assert gw.size % 2 == 0  # even
    nw = gw.size // 2
    nt = 2 * nw

    iw = _matsubara_freq_fermion(beta, nw)
    assert iw.shape == (2*nw,)

    # Subtract 1/iw
    gw_subtract = gw - a / iw

    # Change order to positive-w, negative-w
    gw_subtract = numpy.roll(gw_subtract, nw)

    # fermion to full
    gw_full = numpy.zeros(4*nw, dtype=numpy.complex128)
    gw_full[1::2] = gw_subtract[:]

    # FFT
    gt_full = fft.fft(gw_full) / beta
    assert gt_full.shape == (2*nt,)
    assert numpy.all(numpy.abs(gt_full.imag) < 1e-8)  # real
    gt_full = gt_full.real

    # Extract [0:beta], Add -1/2
    gt_fermion = numpy.zeros(nt+1, dtype=numpy.float64)
    gt_fermion[0:nt] = gt_full[0:nt] - a / 2
    gt_fermion[-1] = - a - gt_fermion[0]

    assert gt_fermion.shape == (nt+1,)
    return gt_fermion


def _fft_fermion_t2w(gt, beta):
    """FFT from G(tau) to G(iw)

    Args:
        gt (numpy.ndarray(nt+1)): G(tau), nt=2*nw
        beta (float): Inverse temperature

    Returns:
        numpy.ndarray(2*nw): G(iw) including w>0 and w<0
    """
    assert gt.size % 2 == 1  # odd
    a =  - gt[0] - gt[-1]
    nt = gt.size - 1
    nw = nt // 2

    iw = _matsubara_freq_fermion(beta, nw)
    assert iw.shape == (2*nw,)

    # Subtract -1/2
    gt_subtract = gt[:-1] + a / 2

    # beta to 2*beta
    gt_full = numpy.append(gt_subtract, -gt_subtract)

    # FFT
    gw_full = fft.ifft(gt_full) * beta
    assert gw_full.shape == (4*nw,)

    # Extract fermion
    gw_fermion = gw_full[1::2]
    assert gw_fermion.shape == (2*nw,)

    # Change order to negative-w, positive-w
    gw_fermion = numpy.roll(gw_fermion, nw)

    # Add 1/iw
    gw_fermion += a / iw

    return gw_fermion


def bgf_fourier_w2t(bgf, tail=None):
    """Fourier transform BlockGf from w to t

    Args:
        bgf (BlockGf(GfImFreq)): Block Green's function in imaginary frequency.
        tail (dict(numpy.ndarray), optional): Coefficient matrix for 1/iw tail. Defaults to None.

    Returns:
        BlockGf(GfImTime): Block Green's function in imaginary time.
    """
    assert isinstance(bgf, BlockGf)
    assert isinstance(bgf.mesh, MeshImFreq)
    assert bgf.mesh.statistic == 'Fermion'
    assert bgf.mesh.positive_only() is False

    beta = bgf.mesh.beta

    nw_pm = bgf.mesh.last_index() - bgf.mesh.first_index() + 1
    assert nw_pm % 2 == 0  # even
    nw = nw_pm // 2  # number of w>0
    nt = nw * 2 + 1

    # Set tail
    if tail is None:
        tail = {}
        for name, gf in bgf:
            _, orb1, orb2 = gf.data.shape
            tail[name] = numpy.eye(orb1, orb2, dtype=numpy.float64)
    else:
        assert isinstance(tail, dict)
        for name, gf in bgf:
            _, orb1, orb2 = gf.data.shape
            assert name in tail
            assert isinstance(tail[name], numpy.ndarray)
            assert tail[name].shape == (orb1, orb2)

    # Make BlockGf
    gf_struct = {name: gf.indices for name, gf in bgf}
    bgf_t = make_block_gf(GfImTime, gf_struct, beta, nt)

    # FFT
    for name, gf in bgf:
        nw_2, norb1, norb2 = gf.data.shape
        assert nw_2 == nw_pm
        assert bgf_t[name].data.shape == (nt, norb1, norb2)
        for i, j in product(range(norb1), range(norb2)):
            gt = _fft_fermion_w2t(gf.data[:, i, j], beta, a=tail[name][i, j])
            assert gt.shape == (nt,)
            bgf_t[name].data[:, i, j] = gt

    return bgf_t
