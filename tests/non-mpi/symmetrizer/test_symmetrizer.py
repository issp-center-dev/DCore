import numpy
from dcore.symmetrizer import pm_symmetrizer
from dcore._dispatcher import BlockGf, GfImFreq
from dcore.tools import make_block_gf, pauli_matrix
import pytest

def _mk_rnd_blockgf(bnames, block_size):
    beta = 1.0
    n_points = 5
    gf_struct = {bname: numpy.arange(block_size) for bname in bnames}
    bgf = make_block_gf(GfImFreq, gf_struct, beta, n_points)
    for _, gf in bgf:
        gf.data[...] = numpy.random.randn(*gf.data.shape)
    return bgf


@pytest.mark.parametrize("transverse", [True, False])
def test_pm_ud(transverse):
    norb = 2
    bgf = _mk_rnd_blockgf(["ud"], 2*norb)
    symm = pm_symmetrizer(norb, True, transverse)
    bgf_symm = symm(bgf)
    orb_I = numpy.identity(norb)
    pauli_mat_ = pauli_matrix()
    if not transverse:
        pauli_mat_ = pauli_mat_[-1:]
    for s_ in pauli_mat_:
        obs = numpy.einsum('st,pq->sptq', s_, orb_I).reshape(2*norb, 2*norb)
        mag = numpy.einsum('wij,ij->w', bgf_symm["ud"].data, obs)
        assert numpy.abs(mag).max() < 1e-10


@pytest.mark.parametrize("transverse", [True, False])
def test_pm_up_down(transverse):
    norb = 2
    bgf = _mk_rnd_blockgf(["up", "down"], norb)
    symm = pm_symmetrizer(norb, False, transverse)
    bgf_symm = symm(bgf)
    orb_I = numpy.identity(norb)
    assert numpy.abs(bgf_symm["up"].data - bgf_symm["down"].data).max() < 1e-10