import irbasis3
import numpy

from dcore.backend.triqs_compat.gf import GfImFreq

default_lambda = 1e+3 # Default value of lambda used for sparse sampling

class Fourier:
    """ Wrapper of """
    pass

def fourier(g, basis=None):
    """ Fourier transform between imaginary-time/frequency domains

    Args:
        g (Gf): Green's function to be transformed
        basis (FiniteTempBasis):
            IR basis
    
    Return: 
        Lazy-evaluation object
    """
    pass


def delta(G0, H0=None, basis=None):
    """ Compute Delta from G0

    Delta(iv) = iv - H0 - G0^{-1}(iv)
    If H0 is None, H0 is obtained by fitting the tail of H0
    """
    assert isinstance(G0, GfImFreq)
    pass


def moments(G, n_moments, basis=None):
    """ Compute the moments for high-frequency expansion
    of fermonic Green's function

    G(iv) = G_1/iv + G_2/(iv)**2 + ..
    by fitting G(iv) with the IR basis

    n_moments: int
        Number of moments to be computed
    return: list of 2D array
        Computed moments [G_1, G_2, ...]
    """
    pass


def dyson(Sigma_iw, G_iw, basis=None):
    pass

def _solve_dyson_for_G0(Sigma_iw, G_iw, basis=None):
    pass


def _solve_dyson_for_Sigma(G0_iw, G_iw, basis=None):
    pass