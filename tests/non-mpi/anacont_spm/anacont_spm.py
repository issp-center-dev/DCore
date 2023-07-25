#
# DCore -- Integrated DMFT software for correlated electrons
# Copyright (C) 2017 The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import numpy as np
from scipy.stats import norm

def _get_matsubara_frequencies_fermionic(n_matsubara, beta):
    return np.pi / beta * (2 * np.arange(0, n_matsubara, 1) + 1)

def _calc_gf_matsubara(energies_dos, dos, beta, n_matsubara):
    matsubara_freq = _get_matsubara_frequencies_fermionic(n_matsubara, beta)
    integrate = lambda wn : np.trapz(y=dos / (1.0j * wn - energies_dos), x=energies_dos)
    gf_matsubara = np.vectorize(integrate)(matsubara_freq)
    return matsubara_freq, gf_matsubara

def _get_dos_semicircular(energies, r):
    return 2.0 / np.pi / r * np.sqrt(np.maximum(1 - np.square(energies / r), 0))

def test_find_sum_rule_const():
    from dcore.anacont_spm import _find_sum_rule_const
    a = 1.123
    b = 2.234
    beta = 40
    n_matsubara = 1000
    n_matsubara_tail = 100
    wn = _get_matsubara_frequencies_fermionic(n_matsubara, beta)
    gf_wn = a - 1j * b / wn
    a_test, b_test = _find_sum_rule_const(wn, gf_wn, n_matsubara_tail, False)
    assert np.allclose(a_test, a, atol=1e-9)
    assert np.allclose(b_test, b, atol=1e-9)

def test_calc_gf_tau_trivial():
    from dcore.anacont_spm import _calc_gf_tau
    a = 1.123
    b = 2.234
    beta = 40
    n_matsubara = 1000
    ntau = 4
    wn = _get_matsubara_frequencies_fermionic(n_matsubara, beta)
    gf_wn = gf_wn = a - 1j * b / wn
    tau_grid, gf_tau = _calc_gf_tau(wn, gf_wn, beta, a, b, ntau)
    tau_grid_expected = np.linspace(0, beta, num=ntau)
    assert np.allclose(tau_grid, tau_grid_expected, atol=1e-10)
    gf_tau_expected = 0.5 * b * np.ones(ntau, dtype=np.float64)
    assert np.allclose(gf_tau, gf_tau_expected, atol=1e-10)

def test_calc_gf_tau_nontrivial():
    from dcore.anacont_spm import _calc_gf_tau, _find_sum_rule_const
    beta = 40
    n_matsubara = 1000
    n_matsubara_tail = 100

    energies_real = np.linspace(-5, 5, num=1000)
    dos = _get_dos_semicircular(energies_real, 4)
    dos_integrated = np.trapz(dos, energies_real)
    assert np.allclose(dos_integrated, 1, atol=1e-13)

    ntau = 11
    wn = _get_matsubara_frequencies_fermionic(n_matsubara, beta)
    _, gf_wn = _calc_gf_matsubara(energies_real, dos, beta, n_matsubara)
    a = 1.123
    gf_wn += a

    a_test, b_test = _find_sum_rule_const(wn, gf_wn, n_matsubara_tail, False)
    assert np.allclose(a_test, a, atol=1e-4)
    assert np.allclose(b_test, 1, atol=1e-3)

    tau_grid, gf_tau = _calc_gf_tau(wn, gf_wn, beta, a_test, b_test, ntau)
    tau_grid_expected = np.linspace(0, beta, num=ntau)
    assert np.allclose(tau_grid, tau_grid_expected, atol=1e-10)
    gf_tau_expected = np.array([0.49991359, 0.04029338, 0.02124658, 0.01544472, 0.0131402, 0.01249759, 0.0131402, 0.01544472, 0.02124658, 0.04029338, 0.49991359])
    assert np.allclose(gf_tau, gf_tau_expected, atol=1e-10)

def test_calc_gf_tau_from_gf_matsubara():
    from dcore.anacont_spm import calc_gf_tau_from_gf_matsubara
    ntau = 11
    n_matsubara = 1000
    n_matsubara_tail = 100
    beta = 40

    energies_real = np.linspace(-5, 5, num=1000)
    dos = _get_dos_semicircular(energies_real, 4)
    dos_integrated = np.trapz(dos, energies_real)
    assert np.allclose(dos_integrated, 1, atol=1e-13)

    ntau = 11
    matsubara_freq, gf_wn = _calc_gf_matsubara(energies_real, dos, beta, n_matsubara)

    tau_grid, gf_tau, const_real_tail, const_imag_tail = calc_gf_tau_from_gf_matsubara(matsubara_freq, gf_wn, ntau, n_matsubara_tail, beta, False)

    tau_grid_expected = np.linspace(0, beta, num=ntau)
    assert np.allclose(const_real_tail, 0, atol=1e-10)
    assert np.allclose(const_imag_tail, 1.0, atol=1e-3)
    assert np.allclose(tau_grid, tau_grid_expected, atol=1e-10)
    gf_tau_expected = np.array([0.49991359, 0.04029338, 0.02124658, 0.01544472, 0.0131402, 0.01249759, 0.0131402, 0.01544472, 0.02124658, 0.04029338, 0.49991359])
    assert np.allclose(gf_tau, gf_tau_expected, atol=1e-10)

def test_get_kernel_matrix():
    from dcore.anacont_spm import _get_kernel_matrix
    beta = 40
    energies = np.linspace(-3, 3, num=5)
    delta_energy = energies[1] - energies[0]
    tau_grid = np.linspace(0, beta, num=3)
    kernel = _get_kernel_matrix(energies, tau_grid, beta, delta_energy)
    kernel_expected = np.array([[5.75073606e-53, 1.31347661e-26, 7.50000000e-01, 1.50000000e+00, 7.50000000e-01], [6.56738307e-27, 1.40364345e-13, 7.50000000e-01, 1.40364345e-13, 6.56738307e-27], [7.50000000e-01, 1.50000000e+00, 7.50000000e-01, 1.31347661e-26, 5.75073606e-53]])
    assert np.allclose(kernel, kernel_expected, atol=1e-10)

def test_get_kernel_matrix_extreme_energies():
    from dcore.anacont_spm import _get_kernel_matrix
    beta = 40
    energies = np.linspace(-30, 30, num=5)
    delta_energy = energies[1] - energies[0]
    tau_grid = np.linspace(0, beta, num=3)
    kernel = _get_kernel_matrix(energies, tau_grid, beta, delta_energy)
    kernel_expected = np.array([
        [0.00000000e+000, 3.97559483e-260, 7.50000000e+000, 1.50000000e+001, 7.50000000e+000], 
        [1.98779741e-260, 7.72230033e-130, 7.50000000e+000, 7.72230033e-130, 1.98779741e-260], 
        [7.50000000e+000, 1.50000000e+001, 7.50000000e+000, 3.97559483e-260, 0.00000000e+000]])
    assert np.allclose(kernel, kernel_expected, atol=1e-10)

def test_getSVD():
    from dcore.anacont_spm import _getSVD, _get_kernel_matrix
    beta = 40
    nsv = 2
    energies = np.linspace(-3, 3, num=5)
    delta_energy = energies[1] - energies[0]
    tau_grid = np.linspace(0, beta, num=3)
    kernel = _get_kernel_matrix(energies, tau_grid, beta, delta_energy)
    U, S, Vt = _getSVD(kernel, nsv=nsv)
    U_expected = np.array([[6.90024281e-01, -7.07106781e-01], [2.18478793e-01, 1.83989109e-16], [6.90024281e-01, 7.07106781e-01]])
    S_expected = np.array([2.02869452, 1.67705098])
    Vt_expected = np.array([[2.55099132e-01, 5.10198264e-01, 5.90968974e-01, 5.10198264e-01, 2.55099132e-01], [3.16227766e-01, 6.32455532e-01, 1.65502277e-16, -6.32455532e-01, -3.16227766e-01]])
    assert np.allclose(U, U_expected, atol=1e-7)
    assert np.allclose(S, S_expected, atol=1e-7)
    assert np.allclose(Vt, Vt_expected, atol=1e-7)

def test_get_svd_for_continuation():
    from dcore.anacont_spm import _get_svd_for_continuation
    beta = 40
    nsv = 2
    ntau = 5
    num_energies = 5

    emin = -4
    emax = 6
    tau_grid = np.linspace(0, beta, num=ntau)

    U, S, Vt, delta_energy, energies_extract = _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies)

    energies_extract_expected = np.linspace(emin, emax, num_energies)
    S_expected = np.array([3.75, 2.79508497])
    U_expected = np.array([
        [-1.00000000e+00, -5.55111512e-17],
        [-2.01777466e-05, -1.08420217e-19],
        [-9.16068221e-10, -7.46069873e-14],
        [-4.15778523e-14, -2.44721856e-07],
        [ 7.80569051e-17, -1.00000000e+00]])
    Vt_expected = np.array([
        [ 2.60189684e-17,  5.20294000e-17, -6.66666667e-01, -6.66666666e-01, -3.33333333e-01],
        [-4.47213595e-01, -8.94427191e-01, -5.34711511e-17, -4.96506831e-17, -2.48253415e-17]])
    assert np.allclose(S, S_expected, atol=1e-7)
    case1 = np.allclose(U, U_expected, atol=1e-7) and np.allclose(Vt, Vt_expected, atol=1e-7)
    case2 = np.allclose(-U, -U_expected, atol=1e-7) and np.allclose(-Vt, -Vt_expected, atol=1e-7)
    assert case1 or case2
    assert np.allclose(delta_energy, (emax - emin) / (num_energies - 1), atol=1e-7)
    assert np.allclose(energies_extract, energies_extract_expected, atol=1e-7)

def test_solveProblem():
    from dcore.anacont_spm import _solveProblem, _find_sum_rule_const, _calc_gf_tau, _get_svd_for_continuation
    beta = 40
    nsv = 100
    n_matsubara = 1000
    n_matsubara_tail = 100
    lambd = 1e-5
    ntau = 1000

    emin = -6
    emax = 6
    num_energies = 1000
    energies = np.linspace(emin, emax, num=num_energies)
    dos = 0.4 * norm.pdf(energies, loc=0.0, scale=0.15)
    dos += 0.3 * norm.pdf(energies, loc=-1.0, scale=0.4)
    dos += 0.3 * norm.pdf(energies, loc=+1.0, scale=0.4)

    wn, gf_wn = _calc_gf_matsubara(energies, dos, beta, n_matsubara)
    a_test, b_test = _find_sum_rule_const(wn, gf_wn, n_matsubara_tail, False)
    tau_grid, gf_tau = _calc_gf_tau(wn, gf_wn, beta, a_test, b_test, ntau)

    U, S, Vt, delta_energy, energies_extract = _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies)

    rho_prime, gf_tau_fit, chi2 = _solveProblem(delta_energy, U, S, Vt, gf_tau, b_test, lambd, verbose=False, max_iters=100, solver='ECOS')
    rho = np.dot(Vt.T, rho_prime)

    d_bhatt = np.trapz(y=np.sqrt(np.multiply(rho, dos)), x=energies_extract)
    assert np.allclose(d_bhatt, 1.0, atol=1e-2)
    assert np.allclose(gf_tau_fit, gf_tau, rtol=1e-3)
    assert np.allclose(chi2, 3.4763690885418575e-06, atol=1e-5)

def test_integral_kramers_kronig():
    from dcore.anacont_spm import _integral_kramers_kronig, dos_to_gf_imag
    from scipy.interpolate import interp1d

    energies = np.linspace(-5, 5, num=1000)
    dos = _get_dos_semicircular(energies, 4)
    gf_imag = dos_to_gf_imag(dos)
    ip = interp1d(energies, gf_imag, fill_value=0, assume_sorted=True)

    results = [
        _integral_kramers_kronig(energies, -3, ip, 1e-10),
        _integral_kramers_kronig(energies, -2, ip, 1e-10),
        _integral_kramers_kronig(energies,  0, ip, 1e-10),
        _integral_kramers_kronig(energies,  2, ip, 1e-10),
        _integral_kramers_kronig(energies,  3, ip, 1e-10)
    ]
    expected_results = [-0.7196181351109775, 
                        -0.4185049136590478, 
                        5.551115123125783e-17, 
                        0.41850491365904774, 
                        0.7196181351109775]
    assert np.allclose(results, expected_results, atol=1e-7)

def test_get_kramers_kronig_realpart():
    from dcore.anacont_spm import get_kramers_kronig_realpart, dos_to_gf_imag
    energies = np.linspace(-5, 5, num=10)
    dos = _get_dos_semicircular(energies, 4)
    gf_imag = dos_to_gf_imag(dos)

    energies_result, gf_real_result, gf_imag_result = get_kramers_kronig_realpart(energies, gf_imag)
    assert np.allclose(energies_result, energies, atol=1e-7)
    assert np.allclose(gf_imag_result, gf_imag, atol=1e-7)
    gf_real_expected = [-0.27741992, -0.38657675, -0.35563321, -0.20813572, -0.06925385, 0.06925385, 0.20813572, 0.35563321, 0.38657675, 0.27741992]
    assert np.allclose(gf_real_result, gf_real_expected, atol=1e-7)

def test_get_single_continuation():
    from dcore.anacont_spm import get_single_continuation, calc_gf_tau_from_gf_matsubara, calc_gf_tau_from_gf_matsubara
    beta = 40
    nsv = 24
    emin = -10
    emax = +10
    num_energies = 100
    lambd = 1e-5
    n_matsubara = 1000
    n_matsubara_tail = 100
    ntau = 25

    energies_dos = np.linspace(-5, 5, num=1000)
    dos = _get_dos_semicircular(energies_dos, 4)
    dos_integrated = np.trapz(dos, energies_dos)
    assert np.allclose(dos_integrated, 1, atol=1e-13)

    matsubara_freq, gf_wn = _calc_gf_matsubara(energies_dos, dos, beta, n_matsubara)

    tau_grid, gf_tau, const_real_tail, const_imag_tail = calc_gf_tau_from_gf_matsubara(matsubara_freq, gf_wn, ntau, n_matsubara_tail, beta, False)

    rho, gf_tau_fit, energies_extract, rho_integrated, chi2 = get_single_continuation(tau_grid, gf_tau, nsv, beta, emin, emax, num_energies, const_imag_tail, lambd, verbose=False)

    assert np.allclose(energies_extract, np.linspace(emin, emax, num=num_energies), atol=1e-7)
    gf_tau_expected = [0.49991358,0.09333313,0.04778339,0.03272457,0.02555529,0.02121236,0.01818985,0.01594869,0.01426303,0.01302357,0.01217094,0.01167125,0.01150657,0.01167125,0.01217094,0.01302357,0.01426303,0.01594869,0.01818985,0.02121236,0.0255553,0.03272457,0.04778339,0.09333313,0.4999136]
    assert np.allclose(gf_tau_fit, gf_tau_expected, atol=1e-7)
    rho_expected = [1.07775066e-02, 2.15550441e-02, 2.15550870e-02, 2.15551471e-02, 2.15552316e-02, 2.15553497e-02, 2.15555150e-02, 2.15557466e-02, 2.15560709e-02, 2.15565250e-02, 2.15571609e-02, 2.15580514e-02, 2.15592984e-02, 2.15610444e-02, 2.15634895e-02, 2.15669134e-02, 2.15717080e-02, 2.15784219e-02, 2.15878231e-02, 2.16009876e-02, 2.16194217e-02, 2.16452342e-02, 2.16813780e-02, 2.17319868e-02, 2.18028475e-02, 2.19020594e-02, 2.20409578e-02, 2.22354006e-02, 2.25075661e-02, 2.28884565e-02, 2.34213769e-02, 2.41667567e-02, 2.52088009e-02, 2.66646059e-02, 2.86965353e-02, 3.15287768e-02, 3.54689939e-02, 4.09355828e-02, 4.84896786e-02, 5.88674767e-02, 7.29997831e-02, 9.19856400e-02, 1.16940499e-01, 1.48531250e-01, 1.85753022e-01, 2.22880063e-01, 2.42065567e-01, 1.96807911e-01, 1.14185563e-09, 2.18439045e-01, 2.18439043e-01, 1.05573341e-09, 1.96807909e-01, 2.42065565e-01, 2.22880061e-01, 1.85753020e-01, 1.48531248e-01, 1.16940497e-01, 9.19856382e-02, 7.29997814e-02, 5.88674750e-02, 4.84896769e-02, 4.09355810e-02, 3.54689923e-02, 3.15287750e-02, 2.86965336e-02, 2.66646041e-02, 2.52087992e-02, 2.41667550e-02, 2.34213750e-02, 2.28884548e-02, 2.25075642e-02, 2.22353988e-02, 2.20409560e-02, 2.19020576e-02, 2.18028457e-02, 2.17319849e-02, 2.16813762e-02, 2.16452325e-02, 2.16194200e-02, 2.16009858e-02, 2.15878213e-02, 2.15784200e-02, 2.15717062e-02, 2.15669117e-02, 2.15634878e-02, 2.15610426e-02, 2.15592965e-02, 2.15580496e-02, 2.15571593e-02, 2.15565232e-02, 2.15560691e-02, 2.15557448e-02, 2.15555132e-02, 2.15553478e-02, 2.15552298e-02, 2.15551454e-02, 2.15550851e-02, 2.15550421e-02, 1.07775057e-02]
    assert np.allclose(rho, rho_expected, atol=1e-3)
    assert np.allclose(rho_integrated, const_imag_tail, atol=1e-3)
    assert np.allclose(chi2, 3.4611947346021525e-06, atol=1e-6)


def test_get_multiple_continuations():
    from dcore.anacont_spm import get_multiple_continuations, calc_gf_tau_from_gf_matsubara, calc_gf_tau_from_gf_matsubara
    beta = 40
    nsv = 24
    emin = -10
    emax = +10
    num_energies = 100
    lambdas = np.logspace(-10, 3, num=3)
    n_matsubara = 1000
    n_matsubara_tail = 100
    ntau = 25

    energies_dos = np.linspace(-5, 5, num=1000)
    dos = _get_dos_semicircular(energies_dos, 4)
    dos_integrated = np.trapz(dos, energies_dos)
    assert np.allclose(dos_integrated, 1, atol=1e-13)

    matsubara_freq, gf_wn = _calc_gf_matsubara(energies_dos, dos, beta, n_matsubara)

    tau_grid, gf_tau, const_real_tail, const_imag_tail = calc_gf_tau_from_gf_matsubara(matsubara_freq, gf_wn, ntau, n_matsubara_tail, beta, False)

    energies_extract, continued_densities, chi2_values, rho_values = get_multiple_continuations(tau_grid, gf_tau, nsv, beta, emin, emax, num_energies, const_imag_tail, lambdas, verbose=False)

    assert np.allclose(energies_extract, np.linspace(emin, emax, num=num_energies), atol=1e-7)
    chi2_expected = [2.4602174332077003e-06, 7.553104010527706e-05, 0.006548520565913404]
    assert np.allclose(chi2_values, chi2_expected, atol=1e-7)
    assert np.allclose(rho_values, const_imag_tail * np.ones(len(rho_values)), atol=1e-7)
    assert len(lambdas) == len(continued_densities)
    assert np.all([num_energies == len(cd) for cd in continued_densities])
