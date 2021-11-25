import numpy as np


def evalulator_high_freq_moment(basis, n):
    """Compute a high-frequency moment of Green's function

    Construct an evaluator for G_n:
        G(iv) \simeq \sum_{n=1}^N G_n/(iv)^n, where
           G_n = (-1)^n (G^{n-1}(0^+) - G^{n-1}(0^-)).

    Attributes:
    -----------
     - `basis` : IR Basis instance
     - `n` : n (moment)
    
    Return a 1D numpy.ndarray of size N [P], where N is the size of the basis.
    The desired high-frequency moment is approximated by
        G_n \simeq sum_{l=0}^{N=1} P_l gl,
    where gl is the expansion coefficients in IR.
    """
    beta = basis.beta
    stat_sign = -1 if basis.statistics == "F" else 1
    uprime = basis.u.deriv(n-1)
    return ((-1)**n) * (uprime(0) - stat_sign * uprime(beta))


def high_freq_moment(gl, basis, num_moments, axis=0):
    """Compute high-frequency moments from IR coefficients of a Green's function

    G(iv) \simeq \sum_{n=1}^N G_n/(iv)^n, where
       G_n = (-1)^n (G^{n-1}(0^+) - G^{n-1}(0^-)).

    Attributes:
    -----------
     - `gl` : Expansion coefficients of Green's function in IR. Three-/one dimensional array.
     - `basis` : IR Basis instance
     - `num_moments` : Number of moments to be computed (>=1)
     - `axis` : Axis of gl corresponding to IR 
    
    Return list `[G_1, G_2, ...]`,
    where G_n are the computed high-frequency moments (each of them is a numpy.ndarray instance).
    """
    assert gl.ndim in [1, 3]

    res = []
    for n_ in range(num_moments):
        res.append(np.einsum('l,l...->...', evalulator_high_freq_moment(basis, n_+1), gl, optimize=True))
    return res