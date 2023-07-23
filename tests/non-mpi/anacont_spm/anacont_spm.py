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
    ntau = 3
    num_energies = 4

    emin = -5
    emax = 5
    tau_grid = np.linspace(0, beta, num=ntau)

    U, S, Vt, delta_energy, energies_extract = _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies)

    energies_extract_expected = np.linspace(emin, emax, num_energies)
    U_expected = np.array([[-1.20839106e-01, 9.92672106e-01], [2.44249065e-15, 2.94209102e-15], [9.92672106e-01, 1.20839106e-01]])
    assert np.allclose(U, U_expected, atol=1e-7)
    S_expected = np.array([3.72677996, 3.72677996])
    assert np.allclose(S, S_expected, atol=1e-7)
    Vt_expected = np.array([[0.44393646, 0.88787292, -0.10808178, -0.05404089], [0.05404089, 0.10808178, 0.88787292, 0.44393646]])
    assert np.allclose(Vt, Vt_expected, atol=1e-7)
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

test_find_sum_rule_const()
test_calc_gf_tau_trivial()
test_calc_gf_tau_nontrivial()
test_calc_gf_tau_from_gf_matsubara()
test_get_kernel_matrix()
test_getSVD()
test_get_svd_for_continuation()
test_solveProblem()
