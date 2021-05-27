import os
import numpy
from numpy.testing import assert_array_equal, assert_array_almost_equal

from h5 import HDFArchive

from dcore import irbasis_x_available
from dcore.dcore_pre import dcore_pre
from dcore.dcore import dcore
from dcore.dcore_gk import dcore_gk
from dcore.tools import float_to_complex_array, save_Sigma_iw_sh_txt, gf_block_names
from dcore._testing import mk_hr_square, create_random_self_energy, gk_from_w90, gk_tail, check_self_energy

import pytest

np = 1
if 'DCORE_TEST_NUM_PROC' in os.environ:
    np = int(os.environ['DCORE_TEST_NUM_PROC'])


def _write_ini(file, spin_orbit, nso_inequiv_sh, corr_to_inequiv, beta, n_iw, nk1, nk2, gk_smpl_freqs):
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
kanamori = [{"(0.0, 0.0, 0.0)," * n_inequiv_sh}]
nk0 = {nk1}
nk1 = {nk2}
nk2 = 1
spin_orbit = {spin_orbit}

[system]
beta = {beta}
n_iw = {n_iw}
fix_mu = True
mu = 0.0

[impurity_solver]
name = null

[control]
max_step = 1
sigma_mix = 0.0
initial_self_energy = Sigma.txt

[tool]
gk_smpl_freqs = {gk_smpl_freqs}
'''

    with open(file, 'w') as f:
        print(ini_str, file=f)


test_square_params = [
    # spin_orbit = False
    ([1],    False, 'gk_smpl_freqs.txt', [0]),
    ([1, 1], False, 'gk_smpl_freqs.txt', [0, 1]),
    ([1],    False, 'gk_smpl_freqs.txt', [0, 0]),
    ([1, 2], False, 'gk_smpl_freqs.txt', [0, 1]),
    # spin_orbit = True
    ([1],    True,  'gk_smpl_freqs.txt', [0]),
    ([1, 1], True,  'gk_smpl_freqs.txt', [0, 1]),
    ([1, 2], True,  'gk_smpl_freqs.txt', [0, 1]),
    ([1],    True,  'dense', [0]),
    ([1],    True,  'sparse', [0]),
    ([1],    True,  'gk_smpl_freqs.txt', [0, 0]), # Two equivalent correlated shells
]
@pytest.mark.parametrize("norb_inequiv_sh, spin_orbit, gk_smpl_freqs, corr_to_inequiv",
    test_square_params)
def test_square(norb_inequiv_sh, spin_orbit, gk_smpl_freqs, corr_to_inequiv, request):
    """
    Multi-orbital square lattice Hubbard model
    """
    numpy.random.seed(100)
    norb_inequiv_sh = numpy.array(norb_inequiv_sh)
    corr_to_inequiv = numpy.array(corr_to_inequiv)
    nso_inequiv_sh = 2 * norb_inequiv_sh

    ncorr_shell = corr_to_inequiv.size
    nf = numpy.sum([nso_inequiv_sh[corr_to_inequiv[icrsh]] for icrsh in range(ncorr_shell)])
    n_inequiv_sh = numpy.unique(corr_to_inequiv).size
    assert n_inequiv_sh == nso_inequiv_sh.size

    nspin = 2
    beta = 5.0
    n_iw = 100
    nk1 = nk2 = 4
    t = 1.0
    noise = 1.0
    Lambda_IR = 1e+4 # default value of dcore_gk
    cutoff_IR = 1e-5 # default value of dcore_gk

    if gk_smpl_freqs == 'sparse' and not irbasis_x_available:
        return

    block_names = gf_block_names(spin_orbit)

    # Block size of each inequivalent shell
    if spin_orbit:
        dim_sh = nso_inequiv_sh
    else:
        dim_sh = nso_inequiv_sh//nspin

    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    # Make Wannier90 file
    mk_hr_square(nf, t, 'square')

    # Random self-energy
    Sigma_iw_sh = \
        create_random_self_energy(block_names, dim_sh, beta, n_iw, noise)
    
    # Save self-energy
    save_Sigma_iw_sh_txt('Sigma.txt', Sigma_iw_sh, block_names)

    # Generate square.init
    _write_ini('square.ini', spin_orbit, nso_inequiv_sh, corr_to_inequiv, 
        beta, n_iw, nk1, nk2, gk_smpl_freqs)

    # Run dcore_pre & dcore
    dcore_pre("square.ini")
    dcore("square.ini", np)
    
    # Check the self-energy saved by dcore
    check_self_energy('square', Sigma_iw_sh, block_names)

    # Sampling fermionic frequencies for dcore_gk
    if gk_smpl_freqs == 'dense':
        smpl_freqs = numpy.arange(-n_iw, n_iw)
    elif gk_smpl_freqs == 'sparse':
        import irbasis_x, irbasis
        b_xy = irbasis_x.twopoint.TruncatedBasis(irbasis.load('F', Lambda_IR), cutoff=cutoff_IR)
        smpl_freqs = b_xy.sampling_points_matsubara(b_xy.dim()-1)
    else:
        smpl_freqs = numpy.array([-2*n_iw, -n_iw, -n_iw+1, -1, 0, 1, n_iw, n_iw+1, 2*n_iw])
        with open('gk_smpl_freqs.txt', 'w') as f:
            print(smpl_freqs.size, file=f)
            for idx, w in enumerate(smpl_freqs):
                print(idx, w, file=f)

    # Run dcore_gk
    dcore_gk("square.ini", np)

    # Read from square_gk.h5
    h = HDFArchive('./square_gk.h5', 'r')
    gkw_read = float_to_complex_array(h['data'])
    smpl_freqs_read = h['smpl_freqs']

    # Reference data
    gkw_ref = gk_from_w90('square', beta, nk1, nk2, 1, Sigma_iw_sh, smpl_freqs, corr_to_inequiv, 0.0)
    #gkw_ref = gk_square(t, nf, beta, nk1, nk2, Sigma_iw_sh, smpl_freqs, corr_to_inequiv)

    assert_array_equal(smpl_freqs, smpl_freqs_read)

    # Sampling frequencies within the window
    lookup_win = numpy.abs(smpl_freqs) < n_iw
    gkw_read = gkw_read.reshape(gkw_ref.shape)
    #for idx_w in range(len(smpl_freqs)):
        #if not lookup_win[idx_w]:
            #continue
        ##print(idx_w, gkw_ref[idx_w, ...].shape, gkw_read[idx_w, ...].shape)
        #print(idx_w, gkw_ref[idx_w, ..., 0])
        #print(idx_w, gkw_read[idx_w, ..., 0])
        #assert_array_almost_equal(gkw_ref[idx_w, ..., 0], gkw_read[idx_w, ..., 0])
    assert_array_almost_equal(gkw_ref[lookup_win, ...].ravel(), gkw_read[lookup_win, ...].ravel())

    # Sampling frequencies outside the window
    lookup_tail = numpy.abs(smpl_freqs) >= n_iw
    gkw_tail = gk_tail(nf, beta, nk1*nk2, smpl_freqs)
    assert_array_almost_equal(gkw_tail[lookup_tail, ...].ravel(), gkw_read[lookup_tail, ...].ravel())

    os.chdir(org_dir)
