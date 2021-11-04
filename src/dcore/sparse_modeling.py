import irbasis3
import numpy

def _compute_sve(lambda_, eps):
    assert eps > 2e-8
    return irbasis3.sve.compute(
        irbasis3.KernelFFlat(lambda_ = lambda_),
        eps = eps, work_dtype = numpy.float64)


# Cache of SVE results
# key: (lambda_, eps)
# value: results of SVE
default_params = (1e+3, 1e-5)
sve_results = {
    default_params: _compute_sve(*default_params)
}

def sve(lambda_, eps):
    """ Get results of SVE """
    key = (lambda_, eps)
    if key not in sve_results:
        sve_results[key] = _compute_sve(*key)
    return sve_results[key]


def finite_basis(beta, statistics='F', lambda_=default_params[0], eps=default_params[1]):
    """ Return a FiniteTempBasis object """
    K = irbasis3.KernelFFlat(lambda_ = lambda_)
    return irbasis3.FiniteTempBasis(
        K, statistics=statistics, beta=beta,
        sve_result = sve(lambda_, eps))

def evalator_tau(basis, tau):
    """ Evaluator at tau"""
    assert isinstance(basis, irbasis3.FiniteTempBasis)
    assert isinstance(tau, numpy.ndarray)

class FourierTransform:
    """ Frourier transform for Green's functions """
    def __init__(self, beta, lambda_, statistics, eps=1e-5):
        self._basis = finite_basis(beta, statistics, eps)
    
    #def to_tau
