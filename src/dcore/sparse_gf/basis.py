import irbasis3
import numpy

def _compute_sve(lambda_, eps):
    assert eps > 2e-8
    return irbasis3.sve.compute(
        irbasis3.KernelFFlat(lambda_ = lambda_),
        eps = eps, work_dtype = numpy.float64)

def _hash_array(x):
    """ Compute a hash of a ndarray object """
    return hash(x.data.tobytes())

class Cache:
    """ Cache objects """
    def __init__(self):
        self.sve_results = {}
        self.tau_sampling = {}
        self.matsubara_sampling = {}

_global_cache = Cache()

def sve(lambda_, eps, cache=_global_cache):
    """ Get results of SVE """
    key = (lambda_, eps)
    if key not in cache.sve_results:
        cache.sve_results[key] = _compute_sve(*key)
    return cache.sve_results[key]

def _tuple_key_basis(basis):
    return (basis.beta, basis.statistics, basis.kernel.lambda_, basis.size)

def finite_temp_basis(beta, statistics, lambda_, eps, cache=_global_cache):
    """ Return a FiniteTempBasis object """
    K = irbasis3.KernelFFlat(lambda_ = lambda_)
    return irbasis3.FiniteTempBasis(K, statistics=statistics, beta=beta,
        sve_result = sve(lambda_, eps, cache))


def _sampling(basis, cache_obj, compute_f, sampling_points):
    """ Return sampling object """
    if sampling_points is None:
        key = _tuple_key_basis(basis)
    else:
        key = _tuple_key_basis(basis) + (_hash_array(sampling_points),)

    if key not in cache_obj:
        cache_obj[key] = compute_f(basis, sampling_points)
    return cache_obj[key]


def tau_sampling(basis, sampling_points=None, cache=_global_cache):
    """ Return TauSampling object """
    return _sampling(basis, cache.tau_sampling, irbasis3.TauSampling, sampling_points)


def matsubara_sampling(basis, sampling_points=None, cache=_global_cache):
    """ Return MatsubaraSampling object """
    return _sampling(basis, cache.matsubara_sampling, irbasis3.MatsubaraSampling, sampling_points)