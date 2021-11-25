import numpy as np

from .block_gf import BlockGf
from .gf import GfImFreq, GfImTime, GfIR, iOmega_n, Identity

from dcorelib.sparse_gf.high_freq import high_freq_moment
from dcorelib.sparse_gf.basis import tau_sampling, matsubara_sampling, finite_temp_basis

def fit_hermitian_tail(giw, basis=None):
    """ 
    Fit the tail of a Green function using a least-squares fitting procedure
    imposing the symmetry  G_{ij}(iv) = G_{ji}^* (-iv).
    Reimplementation of the correponding triqs function.

    In this function, we use sparse modeling.
    Error is not estimated.
    We asssume:
       G(iv) = 1/iv + H/(iv)^2 + O((iv)^3).

    In pratice, we fit
       iv G(iv) - 1 = H/iv + O((iv)^2).

    This works only if the error in G(iv) decays faster than 1/(iv)^2.

    Args:
        giw (GfImFreq): Green's function to be fitted.
    """
    assert isinstance(giw, GfImFreq)
    if basis is None:
        basis = finite_temp_basis(giw.mesh.beta, giw.mesh.statistic)
    
    gl = fit_by_IR(iOmega_n * giw - Identity, basis)
    tail = high_freq_moment(gl.data, basis, 1)
    tail = [0.5*(x + x.T.conjugate()) for x in tail]
    nso = gl.data.shape[1]
    return [np.zeros((nso, nso), dtype=np.complex128), np.identity(nso, dtype=np.complex128)] + tail, 0.0

def delta(giw, basis=None):
    if isinstance(giw, BlockGf):
        delta_iw = giw.copy()
        for bl, g in giw:
            delta_iw[bl] = delta(g, basis)
        return delta_iw
    assert isinstance(giw, GfImFreq), f"Invalid type {type(giw)}"
    if basis is None:
        basis = finite_temp_basis(giw.mesh.beta, giw.mesh.statistic)
    return _delta(giw, basis=basis)

def inverse(g):
    """
    Compute inverse of imaginary-frequency Green's function
    """
    return g.inverse()

def _basis(basis, beta, statistics):
    if basis is not None:
        return basis
    return finite_temp_basis(beta, statistics[0])


def fourier(g, basis=None):
    """ Fourier transform between imaginary-time/frequency domains

    Args:
        g (Gf): Green's function to be transformed
        basis (FiniteTempBasis):
            IR basis
    
    Return: 
        GfIR (fitted result by IR)
    """
    assert type(g) in [GfImFreq, GfImTime]
    basis = _basis(basis, g.mesh.beta, g.mesh.statistic)
    if isinstance(g, GfImFreq):
        smpl = matsubara_sampling(basis, sampling_points=g.mesh.points)
    else:
        smpl = tau_sampling(basis, sampling_points=g.mesh.points)
    return GfIR(smpl.fit(g.data, axis=0), basis)

fit_by_IR = fourier
Fourier = fourier


def _delta(G0, H0=None, basis=None):
    """ Compute Delta from G0 
    Solve
        Delta(iv) = iv - H0 - G0^{-1}(iv).
    If H0 is None, H0 is obtained by fitting the tail of H0 using 
        G0 \simqe I/iv + H0/(iv)^2 + ...
    """
    assert isinstance(G0, GfImFreq)
    if H0 is None:
        mom, _ = fit_hermitian_tail(G0, basis)
        H0 = mom[2]

    delta_iw = G0.copy()
    delta_iw << iOmega_n - H0 - inverse(G0)
    return delta_iw


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
    return high_freq_moment(G.data, basis, n_moments, axis=0)


def dyson(Sigma_iw=None, G_iw=None, G0_iw=None):
    """
    G_iw^-1 = G0_iw^-1 - Sigma_iw
    """
    if Sigma_iw is None:
        return G0_iw.inverse() - G_iw.inverse()
    elif G_iw is None:
        return (G0_iw.inverse() - Sigma_iw).inverse()
    elif G0_iw is None:
        return (G_iw.inverse() + Sigma_iw).inverse()
    else:
        raise RuntimeError("Invalid arguments!")
    