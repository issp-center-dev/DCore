import pytest
import numpy as np
from dcorelib.triqs_compat.utility.pade_approximants import PadeApproximant


def test_pade():
    # C(z_i) = u_i i=1, ..., N
    z = np.array([1j, 2j])
    u = np.array([0.2j, 0.1j])
    N = z.size 
    pade = PadeApproximant(z, u)

    u_reconst = np.array([pade(z_) for z_ in z])
    np.testing.assert_allclose(u, u_reconst)

    u_reconst2 = pade(z)
    np.testing.assert_allclose(u, u_reconst2)