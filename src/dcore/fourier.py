import numpy
from scipy import fft



def _matsubara_freq_fermion(beta, nw):
    iw_positive = numpy.array([(2*i+1) * numpy.pi / beta for i in range(nw)])
    iw = numpy.append(-iw_positive[::-1], iw_positive) * 1j
    return iw


def _fft_fermion_w2t(gw, beta, a=1):
    # nw_2 = gw.size
    assert gw.size % 2 == 0  # even
    nw = gw.size // 2
    nt = 2 * nw
    # print("nt =", nt)
    # print("nw =", nw)

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
    # print(gt_full.shape)
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
    assert gt.size % 2 == 1  # odd
    a =  - gt[0] - gt[-1]
    # print("a=", a)
    nt = gt.size - 1
    nw = nt // 2
    # print("nt =", nt)
    # print("nw =", nw)

    iw = _matsubara_freq_fermion(beta, nw)
    assert iw.shape == (2*nw,)

    # Subtract -1/2
    gt_subtract = gt[:-1] + a / 2
    # print(gt_subtract.shape)

    # beta to 2*beta
    gt_full = numpy.append(gt_subtract, -gt_subtract)
    # print(gt_full.shape)
    # print(gt_full)

    # FFT
    # gw_full = fft.fft(gt_full, norm='backward') * beta
    gw_full = fft.ifft(gt_full) * beta
    # gw_full = fft.ifft(gt_full) * (beta / nt)
    # print("gw_full.shape =", gw_full.shape)
    assert gw_full.shape == (4*nw,)

    # Extract fermion
    gw_fermion = gw_full[1::2]
    assert gw_fermion.shape == (2*nw,)

    # Change order to negative-w, positive-w
    gw_fermion = numpy.roll(gw_fermion, nw)

    # Add 1/iw
    gw_fermion += a / iw

    return gw_fermion



def fermion_fourier(gf, iw, beta, coef_inv_iw=1, direction='w2t'):
    isinstance(gf, numpy.ndarray)
    isinstance(iw, numpy.ndarray)
    assert gf.ndim == 1
    assert iw.ndim == 1
    assert gf.size == iw.size

    n = gf.size


def bgf_fourier(bgf, direction='w2t'):
    isinstance(bgf, BlockGf)

    iw = numpy.array([g.mesh(n) for n in range(g.mesh.first_index(), g.mesh.last_index()+1)])
