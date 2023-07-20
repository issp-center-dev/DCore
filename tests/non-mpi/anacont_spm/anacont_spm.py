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
    kernel_expected = np.array([[1.35363854e-35, 3.09173043e-09, 1.76538950e+17, 3.53077900e+17, 1.76538950e+17], [6.06658705e-35, 6.54519338e-09, 1.76538950e+17, 1.66782191e+17, 3.93911642e+16], [2.71885569e-34, 1.38561745e-08, 1.76538950e+17, 7.87823284e+16, 8.78935678e+15]])
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
    U_expected = np.array([[-0.82897594, -0.5160465], [-0.46366012, 0.41851176], [-0.31275899, 0.74735796]])
    S_expected = np.array([5.13899032e+17, 1.43721900e+17])
    Vt_expected = np.array([[-2.42040421e-52, -1.93255052e-26, -5.51498911e-01, -7.67978126e-01, -3.25666315e-01], [1.54186536e-51, 8.00106734e-26, 7.98202003e-01, -3.72425544e-01, -4.73468875e-01]])
    assert np.allclose(U, U_expected, atol=1e-10)
    assert np.allclose(S, S_expected, atol=1e-10)
    assert np.allclose(Vt, Vt_expected, atol=1e-10)

def test_get_svd_for_continuation():
    from dcore.anacont_spm import _get_svd_for_continuation, _find_sum_rule_const, _calc_gf_tau
    beta = 40
    nsv = 2
    ntau = 3
    num_energies = 4

    emin = -5
    emax = 5
    tau_grid = np.linspace(0, beta, num=ntau)

    U, S, Vt, delta_energy, energies_extract = _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies)

    energies_extract_expected = np.linspace(emin, emax, num_energies)
    U_expected = np.array([[0.9273743, 0.37279142], [0.34485335, -0.818888], [0.14509679, -0.43640464]])
    assert np.allclose(U, U_expected, atol=1e-10)
    S_expected = np.array([9.44496948e+17, 1.29378052e+17])
    assert np.allclose(S, S_expected, atol=1e-10)
    Vt_expected = np.array([[1.53262951e-86, 2.30426322e-29, 9.17662575e-01, 3.97360540e-01], [-3.12088352e-85, -2.58298279e-28, -3.97360540e-01, 9.17662575e-01]])
    assert np.allclose(Vt, Vt_expected, atol=1e-10)
    assert np.allclose(delta_energy, (emax - emin) / (num_energies - 1), atol=1e-10)
    assert np.allclose(energies_extract, energies_extract_expected, atol=1e-10)

def test_solveProblem():
    from dcore.anacont_spm import _solveProblem, _getSVD, _get_kernel_matrix, _find_sum_rule_const, _calc_gf_tau
    beta = 40
    nsv = 35
    n_matsubara = 1000
    n_matsubara_tail = 100
    lambd = 1e-7
    ntau = 1000

    energies = np.linspace(-5, 5, num=1000)
    delta_energy = energies[1] - energies[0]
    dos = _get_dos_semicircular(energies, 4)

    wn, gf_wn = _calc_gf_matsubara(energies, dos, beta, n_matsubara)
    a_test, b_test = _find_sum_rule_const(wn, gf_wn, n_matsubara_tail, False)
    tau_grid, gf_tau = _calc_gf_tau(wn, gf_wn, beta, a_test, b_test, ntau)

    kernel = _get_kernel_matrix(energies, tau_grid, beta, delta_energy)
    U, S, Vt = _getSVD(kernel, nsv=nsv)

    #rho_prime, gf_tau_fit, chi2 = _solveProblem(delta_energy, U, S, Vt, gf_tau, b_test, lambd, verbose=True, max_iters=100, solver='ECOS')
    #assert np.allclose(chi2, 1.0, atol=1e-10)

test_find_sum_rule_const()
test_calc_gf_tau_trivial()
test_calc_gf_tau_nontrivial()
test_calc_gf_tau_from_gf_matsubara()
test_get_kernel_matrix()
test_getSVD()
test_get_svd_for_continuation()
test_solveProblem()
