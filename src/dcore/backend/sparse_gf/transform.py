from dcore.backend.triqs_compat.gf import GfImFreq, GfImTime, GfIR, iOmega_n
from dcore.backend.triqs_compat.gf.gf import inverse
from .basis import tau_sampling, matsubara_sampling, finite_temp_basis
from .high_freq import high_freq_moment


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
    basis = _basis(basis, g.beta, g.statistic)
    if isinstance(g, GfImFreq):
        smpl = matsubara_sampling(basis, sampling_points=g.mesh.points)
    else:
        smpl = tau_sampling(basis, sampling_points=g.mesh.points)
    return GfIR(smpl.fit(g.data, axis=0), basis)



def delta(G0, H0=None, basis=None):
    """ Compute Delta from G0 
    Solve
        Delta(iv) = iv - H0 - G0^{-1}(iv).
    If H0 is None, H0 is obtained by fitting the tail of H0 using 
        G0 \simqe I/iv + H0/(iv)^2 + ...
    """
    assert isinstance(G0, GfImFreq)
    if H0 is None:
        _, H0 = moments(G0, 2, basis=basis)

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
    if Sigma_iw is None:
        return G0_iw.inverse() - G_iw.inverse()
    elif G_iw is None:
        return (G0_iw.inverse() - Sigma_iw).inverse()
    else:
        raise RuntimeError("Invalid arguments!")