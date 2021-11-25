from dcorelib.sparse_gf.basis import *
import numpy
import pytest

@pytest.mark.parametrize("statistics", ['F', 'B'])
def test_basis(statistics):
    lambda_ = 100
    beta = 2
    eps = 1e-5
    basis = finite_temp_basis(beta, statistics, lambda_, eps)
    assert basis.statistics == statistics
    assert basis.wmax == lambda_/beta

    gl = numpy.random.randn(basis.size)

    tau_smpl = tau_sampling(basis)
    tau_smpl2 = tau_sampling(basis)
    numpy.testing.assert_array_equal(tau_smpl.evaluate(gl), tau_smpl2.evaluate(gl))

    matsubara_smpl = matsubara_sampling(basis)
    matsubara_smpl2 = matsubara_sampling(basis)
    numpy.testing.assert_array_equal(matsubara_smpl.evaluate(gl), matsubara_smpl2.evaluate(gl))
