import irbasis3
import numpy

from dcore.triqs_compat.gf import GfImFreq


def delta(G0, H0=None):
    """ Compute Delta from G0

    Delta(iv) = iv - H0 - G0^{-1}(iv)
    If H0 is None, H0 is obtained by fitting the tail of H0
    """
    assert isinstance(G0, GfImFreq)
    pass


def inverse_dyson(Sigma_iw, G_iw):
    pass


def dyson(G0_iw, G_iw):
    pass