import numpy
import irbasis
import irbasis_x

def construct_basis(stat, beta, Lambda, cutoff):
    """
    stat: str
        'F' or 'B'
    beta: float
        Inverse temperature
    Lambda: float
        IR Lambda
    cutoff: float
        IR cutoff parameter
    """
    b_xy = irbasis_x.twopoint.TruncatedBasis(irbasis.load(stat, Lambda), cutoff=cutoff)
    return irbasis_x.twopoint.FiniteTemperatureBasis(b_xy, beta)

def fit_iw(basis, gw, axis=0):
    """
    Fit Green's function on sampling frequencies by IR

    basis: FiniteTemperatureBasis
      Two-point basis
    gw: ndarray
      Green's function on sampling frequencies
    axis: int
      Position of the axis for the frequencies
    """
    if gw.ndim == 1:
        return basis.fit_iw(gw)
    gw = numpy.moveaxis(gw, [axis], [0])
    nw = gw.shape[0]
    rest_shapes = gw.shape[1:]
    gl = basis.fit_iw(gw.reshape(nw, -1))
    return numpy.moveaxis(gl.reshape((basis.dim(),) + rest_shapes), [0], [axis])

