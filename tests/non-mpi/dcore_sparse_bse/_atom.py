from irbasis_x import atom
from triqs.gf import GfImFreq, BlockGf
import numpy

def sigma(U, beta, n_iw, spin_orbit):
    vsample = 2 * numpy.arange(-n_iw, n_iw) + 1
    sigma_data = atom.sigma(U, beta, vsample)

    if spin_orbit:
        g = GfImFreq(indices = [0, 1], beta = beta, n_points = n_iw)
        g.data[:, :, :] = 0.0
        g.data[:, 0, 0] = sigma_data
        g.data[:, 1, 1] = sigma_data
        return BlockGf(name_list = ('ud',), block_list = (g,), make_copies = False)
    else:
        g = GfImFreq(indices = [0], beta = beta, n_points = n_iw)
        g.data[:, 0, 0] = sigma_data
        return BlockGf(name_list = ('up', 'down'), block_list = (g,g), make_copies = True)

def _g2loc_uu_ud(U, beta, v, vp, w):
    chi_uu, chi_ud = atom.chi_ph(U, beta, v, vp, w)
    gfat = lambda n: atom.gf(U, beta, n)
    discon = beta * atom._delta(w, 0) * gfat(w + v) * gfat(vp)
    return chi_uu + discon, chi_ud + discon

def G2loc(U, beta, wsample_ph):
    v, vp, w = wsample_ph

    #chi_uu, chi_ud = atom.chi_ph(U, beta, *wsample_ph)
    g2loc_uu, g2loc_ud = _g2loc_uu_ud(U, beta, *wsample_ph)
    up, dn = 0, 1
    g2loc = numpy.zeros((len(wsample_ph[0]), 2, 2, 2, 2), dtype=numpy.complex128)
    #g2loc[:, up, up, up, up] = g2loc[:, dn, dn, dn, dn] = chi_uu + discon
    #g2loc[:, up, up, dn, dn] = g2loc[:, dn, dn, up, up] = chi_ud + discon
    g2loc[:, up, up, up, up] = g2loc[:, dn, dn, dn, dn] = g2loc_uu
    g2loc[:, up, up, dn, dn] = g2loc[:, dn, dn, up, up] = g2loc_ud

    _, g2loc_ud_cr = _g2loc_uu_ud(U, beta, vp+w, vp, v-vp)
    g2loc[:, up, dn, dn, up] = g2loc[:, dn, up, up, dn] = -g2loc_ud_cr

    return g2loc