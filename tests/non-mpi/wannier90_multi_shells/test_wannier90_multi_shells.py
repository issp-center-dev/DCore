import os
import numpy
from numpy.testing import assert_allclose

from dcorelib.triqs_compat.h5 import HDFArchive

from dcore.dcore_pre import dcore_pre
from dcore.dcore import dcore
from dcore.tools import float_to_complex_array, save_Sigma_iw_sh_txt, gf_block_names
from dcore._testing import mk_hr_square, create_random_self_energy

import pytest


def _write_ini(file, spin_orbit, nso_corr_sh, nf, beta, n_iw, nk1, nk2, Lambda, cutoff):
    norb = nf//2
    n_corr_sh = n_inequiv_sh = len(nso_corr_sh)
    n_orb_str = ' '.join(map(str, (nso_corr_sh//2).tolist()))

    ini_str = \
    f'''[model]
lattice = wannier90
seedname = square
t = 1.0
norb = {n_orb_str}
ncor = {n_corr_sh}
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
'''

# [tool]
# Lambda_IR = {Lambda}
# cutoff_IR = {cutoff}'''

    with open(file, 'w') as f:
        print(ini_str, file=f)




test_square_params = [
    (numpy.array([2]),    True),
    (numpy.array([2, 2]), True),
    (numpy.array([2, 4]), True),
    (numpy.array([6, 8]), True),
    (numpy.array([2, 2]), False),
    #(numpy.array([2, 4]),  False),
]
@pytest.mark.parametrize("nso_corr_sh, spin_orbit", test_square_params)
def test_square(nso_corr_sh, spin_orbit, request):
    """
    Each orbital forms a correlated shell
    """
    numpy.random.seed(100)

    assert all(nso_corr_sh % 2 == 0)

    nf = numpy.sum(nso_corr_sh)
    ncorr_shell = n_inequiv_sh = len(nso_corr_sh)

    nspin = 2
    beta = 5.0
    n_iw = 1000
    nk1 = nk2 = 4
    Lambda = 1e+4
    cutoff = 1e-5
    t = 1.0

    # Block size of each inequivalent shell
    if spin_orbit:
        dim_sh = nso_corr_sh
    else:
        dim_sh = nso_corr_sh//nspin
    noise = 0.1

    block_names = gf_block_names(spin_orbit)

    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    # Make Wannier90 file
    if spin_orbit:
        mk_hr_square(nf, t, 'square')
    else:
        mk_hr_square(nf//2, t, 'square')

    # Random self-energy
    Sigma_iw_sh = \
        create_random_self_energy(block_names, dim_sh, beta, n_iw, noise)
    
    # Save self-energy
    for ish in range(len(Sigma_iw_sh)):
        print(ish, Sigma_iw_sh[ish])
    save_Sigma_iw_sh_txt('Sigma.txt', Sigma_iw_sh, block_names)

    # Generate square.init
    _write_ini('square.ini', spin_orbit, nso_corr_sh, nf, beta, n_iw, nk1, nk2, Lambda, cutoff)

    # Run dcore_pre & dcore
    dcore_pre("square.ini")
    dcore("square.ini")
    
    # Check the self-energy saved by dcore
    with HDFArchive('square.out.h5', 'r') as h:
        for ish in range(ncorr_shell):
            for bname in block_names:
                x = float_to_complex_array(h['dmft_out']['Sigma_iw']['ite1']
                    [f"""sh{ish}"""][bname]['data'])
                y = Sigma_iw_sh[ish][bname].data
                assert_allclose(x, y)

    os.chdir(org_dir)
