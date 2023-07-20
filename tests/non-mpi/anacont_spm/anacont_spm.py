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
    matsubara_freq, gf_wn = _calc_gf_matsubara(energies_real, dos, beta, n_matsubara)
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
    wn = _get_matsubara_frequencies_fermionic(n_matsubara, beta)
    matsubara_freq, gf_wn = _calc_gf_matsubara(energies_real, dos, beta, n_matsubara)

    tau_grid, gf_tau, const_real_tail, const_imag_tail = calc_gf_tau_from_gf_matsubara(matsubara_freq, gf_wn, ntau, n_matsubara_tail, beta, False)

    tau_grid_expected = np.linspace(0, beta, num=ntau)
    assert np.allclose(const_real_tail, 0, atol=1e-10)
    assert np.allclose(const_imag_tail, 1.0, atol=1e-3)
    assert np.allclose(tau_grid, tau_grid_expected, atol=1e-10)
    gf_tau_expected = np.array([0.49991359, 0.04029338, 0.02124658, 0.01544472, 0.0131402, 0.01249759, 0.0131402, 0.01544472, 0.02124658, 0.04029338, 0.49991359])
    assert np.allclose(gf_tau, gf_tau_expected, atol=1e-10)

test_find_sum_rule_const()
test_calc_gf_tau_trivial()
test_calc_gf_tau_nontrivial()
test_calc_gf_tau_from_gf_matsubara()
