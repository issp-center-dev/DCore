import numpy
from dcore.symmetrizer import pm_symmetrizer
from dcore._dispatcher import BlockGf, GfImFreq
from dcore.tools import make_block_gf, pauli_matrix

def _mk_rnd_blockgf(bnames, block_size):
    beta = 1.0
    n_points = 5
    gf_struct = {bname: numpy.arange(block_size) for bname in bnames}
    bgf = make_block_gf(GfImFreq, gf_struct, beta, n_points)
    for _, gf in bgf:
        gf.data[...] = numpy.random.randn(*gf.data.shape)
    return bgf


def test_pm_ud():
    norb = 2
    bgf = _mk_rnd_blockgf(["ud"], 2*norb)
    symm = pm_symmetrizer(norb, True)
    bgf_symm = symm(bgf)
    orb_I = numpy.identity(norb)
    for s_ in pauli_matrix():
        obs = numpy.einsum('st,pq->sptq', s_, orb_I).reshape(2*norb, 2*norb)
        mag = numpy.einsum('wij,ij->w', bgf_symm["ud"].data, obs)
        assert numpy.abs(mag).max() < 1e-10

def test_pm_up_down():
    norb = 2
    bgf = _mk_rnd_blockgf(["up", "down"], norb)
    symm = pm_symmetrizer(norb, False)
    bgf_symm = symm(bgf)
    orb_I = numpy.identity(norb)
    assert numpy.abs(bgf_symm["up"].data - bgf_symm["down"].data).max() < 1e-10