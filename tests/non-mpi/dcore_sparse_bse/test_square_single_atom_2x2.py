import os
import numpy
from numpy.testing import assert_allclose, assert_array_almost_equal

from h5 import HDFArchive

from dcore import irbasis_x_available
from dcore.dcore_pre import dcore_pre
from dcore.dcore import dcore
from dcore.dcore_gk import dcore_gk
from dcore.dcore_sparse_bse import dcore_sparse_bse
from dcore.dcore_vertex import dcore_vertex

from dcore.tools import float_to_complex_array, gf_block_names
from dcore._testing import mk_hr_square_2x2, gk_tail, check_self_energy, gk_from_w90

import pytest

from dcore import irbasis_x_available

np = 1
if 'DCORE_TEST_NUM_PROC' in os.environ:
    np = int(os.environ['DCORE_TEST_NUM_PROC'])


def _write_ini(file, spin_orbit, nso_inequiv_sh, corr_to_inequiv, U, beta,
    n_iw, nk1, nk2, gk_smpl_freqs, Lambda_IR, cutoff_IR, cutoff_IR_2P):
    n_corr_sh = corr_to_inequiv.size
    n_inequiv_sh = nso_inequiv_sh.size
    n_orb_str = ' '.join(map(str, (nso_inequiv_sh//2).tolist()))
    corr_to_inequiv_str = ' '.join(map(str, corr_to_inequiv))

    ini_str = \
    f'''[model]
lattice = wannier90
seedname = square
t = 1.0
norb = {n_orb_str}
ncor = {n_corr_sh}
corr_to_inequiv = {corr_to_inequiv_str}
interaction = kanamori
kanamori = [{f"""({U}, 0.0, 0.0),""" * n_inequiv_sh}]
nk0 = {nk1}
nk1 = {nk2}
nk2 = 1
spin_orbit = {spin_orbit}

[system]
beta = {beta}
n_iw = {n_iw}
fix_mu = True
mu = {0.5*U}

[impurity_solver]
name = pomerol
exec_path{{str}}=pomerol2dcore

[control]
max_step = 1
sigma_mix = 1.0
#initial_self_energy = Sigma.txt

[tool]
gk_smpl_freqs = {gk_smpl_freqs}
Lambda_IR = {Lambda_IR}
cutoff_IR = {cutoff_IR}
cutoff_IR_2P = {cutoff_IR_2P}

[vertex]
wb_ph = 0

[sparse_bse]
qsample = qsample.txt
'''

    with open(file, 'w') as f:
        print(ini_str, file=f)


test_square_params = [
    # spin_orbit = False
    ([1],    False, [0, 0, 0, 0]),
    ([1, 1, 1, 1],    False, [0, 1, 2, 3]),
    # spin_orbit = True
    ([1],    True, [0, 0, 0, 0]),
    ([1, 1, 1, 1],    True, [0, 1, 2, 3]),
]
@pytest.mark.skipif(not irbasis_x_available, reason="irbasis_x is not available")
@pytest.mark.parametrize("norb_inequiv_sh, spin_orbit, corr_to_inequiv",
    test_square_params)
def test_square(norb_inequiv_sh, spin_orbit, corr_to_inequiv, request):
    """
    Single-orbital Hubbard model with atomic vertex
    A unit cell may contain multiple correlated shells, but the transfer integral is digonal
    in terms of shells.
    """
    from dcore.irbasis_util import construct_basis
    from _atom import sigma, G2loc

    numpy.random.seed(100)
    norb_inequiv_sh = numpy.array(norb_inequiv_sh)
    corr_to_inequiv = numpy.array(corr_to_inequiv)
    nso_inequiv_sh = 2 * norb_inequiv_sh

    ncorr_shell = corr_to_inequiv.size
    nf = numpy.sum([nso_inequiv_sh[corr_to_inequiv[icrsh]] for icrsh in range(ncorr_shell)])
    n_inequiv_sh = numpy.unique(corr_to_inequiv).size
    assert n_inequiv_sh == nso_inequiv_sh.size

    W = 8.0 # Band width
    U = 2.5 * W # Onsite repulsion
    mu = 0.5 * U
    TN = 4/U # Strong-coupling limit of TN (4t^2/U)
    beta = 1.05/TN # Just below TN
    n_iw = 1000
    nk1 = nk2 = 5
    nspin = 2 if spin_orbit else 1
    t = 1.0
    Lambda_IR = 1e+2 # default value of dcore_gk
    cutoff_IR = 1e-5 # default value of dcore_gk
    assert Lambda_IR/beta > W

    block_names = gf_block_names(spin_orbit)

    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    # Run dcore_pre
    mk_hr_square_2x2(nspin, t, 'square')
    _write_ini('square.ini', spin_orbit, nso_inequiv_sh, corr_to_inequiv, 
        U, beta, n_iw, nk1, nk2, 'sparse', Lambda_IR, cutoff_IR, cutoff_IR)
    dcore_pre("square.ini")

    # Run dcore
    dcore("square.ini")
    Sigma_iw_ref = [sigma(U, beta, n_iw, spin_orbit)] * n_inequiv_sh
    check_self_energy('square', Sigma_iw_ref, block_names)

    # Run dcore_gk
    dcore_gk("square.ini", np)
    basis_f = construct_basis('F', beta, Lambda_IR, cutoff_IR)
    smpl_freqs = basis_f.wsample//2
    with HDFArchive('./square_gk.h5', 'r') as h:
       gkw_read = float_to_complex_array(h['data'])
       #smpl_freqs_read = h['smpl_freqs']
    gkw_ref = gk_from_w90('square', beta, nk1, nk2, 1, Sigma_iw_ref, smpl_freqs, corr_to_inequiv, mu)
    gkw_read = gkw_read.reshape(gkw_ref.shape)

    # Sampling frequencies within the window
    lookup_win = numpy.abs(smpl_freqs) < n_iw
    #print(gkw_ref.shape)
    #print(gkw_read.shape)
    #print(gkw_ref[10, 0, :, 0])
    #print(gkw_read[10, 0, :, 0])
    #print(gkw_read[lookup_win[0], 0, :, 0].shape)
    assert_array_almost_equal(gkw_ref[lookup_win, ...].ravel(), gkw_read[lookup_win, ...].ravel())

    # Sampling frequencies outside the window
    lookup_tail = numpy.abs(smpl_freqs) >= n_iw
    gkw_tail = gk_tail(nf, beta, nk1*nk2, smpl_freqs)
    assert_array_almost_equal(gkw_tail[lookup_tail, ...].ravel(), gkw_read[lookup_tail, ...].ravel())

    # bosonic momenta
    with open('qsample.txt', 'w') as f:
        print(1, file=f)
        print(0, 0, 0, 0, file=f)
    dcore_vertex("square.ini")

    with HDFArchive('square_vertex.h5', 'r') as h:
        wsample_ph_g2loc = h['wsample_ph']
        G2_loc_read = float_to_complex_array(h['G2loc/sh0'])
        G2_loc_ref = G2loc(U, beta, wsample_ph_g2loc)
        assert_allclose(G2_loc_ref, G2_loc_read.reshape(G2_loc_ref.shape))

    dcore_sparse_bse("square.ini", np)

    with HDFArchive('square_chi.h5', 'r') as h:
        # chi: (nwb, nq, nf, nf, nf, nf)
        chi = float_to_complex_array(h['chi'])
    chi_mat = chi[0, 0, :, :, :, :].\
        reshape((nf,)*4).\
        transpose((0,1,3,2)).reshape(nf**2, nf**2)
    atol = numpy.abs(chi_mat).max()
    print(chi_mat)
    print(chi_mat.T.conj())
    assert_allclose(chi_mat, chi_mat.T.conj(), rtol=0, atol=1e-3*atol)
    evals, _ = numpy.linalg.eigh(chi_mat)

    print("eval", evals)
    #assert abs(evals[-1]) < 1e-3  # One irrelevant one
    assert all(evals[:3] < -5e+1) # Correspond to (pi, pi)

    os.chdir(org_dir)