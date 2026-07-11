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

import warnings

import numpy
from scipy import fft
from itertools import product
from dcore._dispatcher import BlockGf, Gf, GfImFreq, GfImTime, MeshImFreq
from dcore.tools import make_block_gf
from dcore import ir_basis


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


def _ir_default_wmax(beta, nw):
    """Fallback real-frequency cutoff = largest Matsubara frequency on the grid.

    This over-estimate guarantees the basis spans everything the sampling grid
    resolves, but it scales with the *number* of frequencies (a numerical
    parameter), not the physical spectral width. For large nw it makes
    Lambda = beta*wmax large, which is slower to build and more ill-conditioned
    (hence less accurate) than a physically-sized cutoff. Emits a warning so the
    fallback is never silent; callers should pass an explicit wmax covering the
    spectral support whenever it is known.

    NOTE for the future consumer wiring: the physically-correct cutoff is the
    dispersion spectral half-range max|eps_k - mu|, which this pure G(iw)<->G(tau)
    transform cannot see (it has no eps(k)/mu). The consumer that owns the lattice
    (DMFT driver / SumkDFT) must supply wmax = max|eps_k - mu| via [system] ir_wmax
    -- do NOT re-derive a band/grid heuristic here: the sister project H-wave hit
    exactly that bug (issp-center-dev/H-wave issue #57), where a naive band measure
    produced an ill-conditioned basis and wrong results.
    """
    wmax = (2 * nw - 1) * numpy.pi / beta
    warnings.warn(
        f"IR basis wmax not specified; falling back to the Matsubara-grid-edge "
        f"wmax={wmax:.3g} (Lambda=beta*wmax={beta * wmax:.3g}), which is slower "
        f"and less accurate than a physically-sized cutoff. Pass an explicit "
        f"wmax (or [system] ir_wmax) covering the spectral support.",
        stacklevel=3,
    )
    return wmax


def _ir_fermion_w2t(gw, beta, wmax=None, eps=1e-10):
    """FFT from G(iw) to G(tau) via the IR (sparse-ir) basis.

    Args:
        gw (numpy.ndarray(2*nw)): G(iw) including w>0 and w<0 on the dense symmetric grid.
        beta (float): Inverse temperature.
        wmax (float, optional): Real-frequency cutoff. Defaults to the largest Matsubara
            frequency on the grid.
        eps (float, optional): Basis truncation tolerance. Defaults to 1e-10.

    Returns:
        numpy.ndarray(nt+1): real G(tau) on linspace(0, beta, nt+1), nt=2*nw.
    """
    import sparse_ir
    assert gw.size % 2 == 0  # even
    nw = gw.size // 2
    nt = 2 * nw
    if wmax is None:
        wmax = _ir_default_wmax(beta, nw)

    basis = ir_basis.get_basis(beta, wmax, eps, 'F')
    # Fermionic Matsubara indices 2n+1 for n in [-nw, nw), matching _matsubara_freq_fermion.
    n_idx = numpy.array([2 * n + 1 for n in range(-nw, nw)])
    smpl_w = sparse_ir.MatsubaraSampling(basis, sampling_points=n_idx)
    tau_grid = numpy.linspace(0.0, beta, nt + 1)
    smpl_t = sparse_ir.TauSampling(basis, sampling_points=tau_grid)

    g_l = smpl_w.fit(gw)
    gt = smpl_t.evaluate(g_l)
    # G(tau) of a fermionic G(iw) is real; mirror the FFT path's loud check
    # rather than silently discarding an imaginary part (a matrix-valued
    # off-diagonal block that is genuinely complex must not pass silently).
    assert numpy.all(numpy.abs(gt.imag) < 1e-8), \
        "IR w2t produced a non-real G(tau); input may be a complex off-diagonal block"
    return gt.real


def _ir_fermion_t2w(gt, beta, wmax=None, eps=1e-10):
    """FFT from G(tau) to G(iw) via the IR (sparse-ir) basis.

    Args:
        gt (numpy.ndarray(nt+1)): real G(tau) on linspace(0, beta, nt+1), nt=2*nw.
        beta (float): Inverse temperature.
        wmax (float, optional): Real-frequency cutoff. Defaults to the largest Matsubara
            frequency on the grid.
        eps (float, optional): Basis truncation tolerance. Defaults to 1e-10.

    Returns:
        numpy.ndarray(2*nw): complex G(iw) including w>0 and w<0 on the dense symmetric grid.
    """
    import sparse_ir
    assert gt.size % 2 == 1  # odd
    nt = gt.size - 1
    nw = nt // 2
    if wmax is None:
        wmax = _ir_default_wmax(beta, nw)

    basis = ir_basis.get_basis(beta, wmax, eps, 'F')
    tau_grid = numpy.linspace(0.0, beta, nt + 1)
    smpl_t = sparse_ir.TauSampling(basis, sampling_points=tau_grid)
    n_idx = numpy.array([2 * n + 1 for n in range(-nw, nw)])
    smpl_w = sparse_ir.MatsubaraSampling(basis, sampling_points=n_idx)

    g_l = smpl_t.fit(gt)
    return smpl_w.evaluate(g_l)


def bgf_fourier_w2t(bgf, tail=None, method='fft', ir_params=None):
    """Fourier transform BlockGf from w to t

    Args:
        bgf (BlockGf(GfImFreq)): Block Green's function in imaginary frequency.
        tail (dict(numpy.ndarray), optional): Coefficient matrix for 1/iw tail. Defaults to None.
        method (str, optional): 'fft' (dense FFT, default) or 'ir' (sparse-ir basis).
        ir_params (dict, optional): {'wmax': float|None, 'eps': float} passed to the IR
            path. Ignored when method='fft'. `tail` is ignored when method='ir'.

    Returns:
        BlockGf(GfImTime): Block Green's function in imaginary time.
    """
    assert isinstance(bgf, BlockGf)
    assert isinstance(bgf.mesh, MeshImFreq)
    assert bgf.mesh.statistic == 'Fermion'
    assert bgf.mesh.positive_only() is False

    if method not in ('fft', 'ir'):
        raise ValueError(f"Unknown method '{method}'; expected 'fft' or 'ir'.")

    beta = bgf.mesh.beta

    nw_pm = bgf.mesh.size
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
            if method == 'fft':
                gt = _fft_fermion_w2t(gf.data[:, i, j], beta, a=tail[name][i, j])
            else:  # method == 'ir'
                gt = _ir_fermion_w2t(gf.data[:, i, j], beta, **(ir_params or {}))
            assert gt.shape == (nt,)
            bgf_t[name].data[:, i, j] = gt

    return bgf_t
