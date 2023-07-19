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
import cvxpy as cp
from scipy.interpolate import interp1d
from scipy.sparse.linalg import svds

def _find_sum_rule_const(matsubara_frequencies, gf_wn, ntail, show_fit=False):
    wntail = matsubara_frequencies[-ntail:]
    gftail = gf_wn[-ntail:]
    const_imagpart = -np.mean(np.imag(gftail) * wntail)
    const_realpart = np.mean(np.real(gftail))
    if show_fit:
        import matplotlib.pyplot as plt
        z = np.linspace(wntail[0], wntail[-1], num=1000)
        plt.scatter(wntail, np.imag(gftail), zorder=5, color='C0', label='data')
        plt.plot(z, -const_imagpart / z, zorder=10, color='C1', label='fit')
        plt.xlabel(r'$\omega_n$')
        plt.ylabel(r'Im $G( \omega_n )$')
        plt.legend()
        plt.tight_layout()
        plt.show()
    return const_realpart, const_imagpart

def _calc_gf_tau(matsubara_frequencies, gf_wn, beta, const_real_tail, const_imag_tail, ntau):
    tau_grid = np.linspace(0, beta, num=ntau)
    gf_tau = -0.5 * const_imag_tail * np.ones(tau_grid.shape, dtype=np.float64)
    kernel = lambda tau : np.sum((np.real(gf_wn) - const_real_tail) * np.cos(matsubara_frequencies * tau) + (np.imag(gf_wn) + const_imag_tail / matsubara_frequencies) * np.sin(matsubara_frequencies * tau))
    gf_tau += 2 / beta * np.vectorize(kernel)(tau_grid)
    gf_tau *= -1
    return tau_grid, gf_tau

def calc_gf_tau_from_gf_matsubara(matsubara_frequencies, gf_wn, ntau, ntail, beta, show_fit=False):
    const_real_tail, const_imag_tail = _find_sum_rule_const(matsubara_frequencies, gf_wn, ntail, show_fit=show_fit)
    print('Determined tail constants: {}, {}'.format(const_real_tail, const_imag_tail))
    tau_grid, gf_tau = _calc_gf_tau(matsubara_frequencies, gf_wn, beta, const_real_tail, const_imag_tail, ntau)
    return tau_grid, gf_tau, const_real_tail, const_imag_tail

def get_single_continuation(tau_grid, gf_tau, nsv, beta, emin, emax, num_energies, sum_rule_const, lambd, verbose=True, max_iters=100, solver='ECOS'):
    U, S, Vt, delta_energy, energies_extract = _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies)
    rho_prime, gf_tau_fit, chi2 = _solveProblem(delta_energy, U, S, Vt, gf_tau, sum_rule_const, lambd, verbose=verbose, max_iters=max_iters, solver=solver)
    rho = np.dot(Vt.T, rho_prime)
    rho_integrated = np.trapz(y=rho, x=energies_extract)
    return rho, gf_tau_fit, energies_extract, rho_integrated, chi2

def get_multiple_continuations(tau_grid, gf_tau, nsv, beta, emin, emax, num_energies, sum_rule_const, lambdas, verbose=True, max_iters=100, solver='ECOS'):
    U, S, Vt, delta_energy, energies_extract = _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies)
    continued_densities = []
    chi2_values = []
    rho_values = []
    for lambd in lambdas:
        rho_prime, _, chi2 = _solveProblem(delta_energy, U, S, Vt, gf_tau, sum_rule_const, lambd, verbose=verbose, max_iters=max_iters, solver=solver)
        rho = np.dot(Vt.T, rho_prime)
        rho_integrated = np.trapz(y=rho, x=energies_extract)
        continued_densities.append(rho)
        chi2_values.append(chi2)
        rho_values.append(rho_integrated)
    return energies_extract, continued_densities, chi2_values, rho_values

def _get_kernel_matrix(energies, tau_grid, beta, delta_energy):
    assert tau_grid[0] == 0
    assert tau_grid[-1] == beta
    den = 1.0 / (1.0 + np.exp(-beta * energies))
    kernel = den * np.exp(np.multiply.outer(-tau_grid, energies))
    kernel[:, 0] *= 0.5
    kernel[:, -1] *= 0.5
    kernel *= delta_energy
    return kernel

def _getSVD(matrix, nsv=None):
    if nsv == None or nsv >= min(matrix.shape) or nsv <= 0:
        raise Exception("nsv must be 0 < nsv < {} (min(matrix.shape))".format(min(matrix.shape)))
    U, S, V = svds(matrix, k=nsv, which='LM')
    S = np.flip(S, axis=0)
    U = np.flip(U, axis=1)
    Vt = np.flip(V, axis=0)
    return U, S, Vt

def _solveProblem(delta_energy, U, S, Vt, gf_tau, sum_rule_const, lambd, verbose=True, max_iters=100, solver='ECOS'):
    Q = len(S)
    Smat = np.diag(S)
    rho_prime = cp.Variable(Q)
    errTerm = gf_tau - U @ Smat @ rho_prime
    objective = cp.Minimize(0.5 * cp.norm2(errTerm) ** 2 + lambd * cp.norm1(rho_prime))
    V_mod = np.copy(Vt.T)
    V_mod[0, :] *= 0.5
    V_mod[-1, :] *= 0.5
    constraints = [Vt.T @ rho_prime >= 0, cp.sum(delta_energy * V_mod @ rho_prime) == sum_rule_const] #uniform real energy grid is assumed here
    prob = cp.Problem(objective, constraints)
    _ = prob.solve(verbose=verbose, solver=solver, max_iters=max_iters)
    gf_tau_fit = np.dot(U, np.dot(Smat, rho_prime.value))
    chi2 = 0.5 * np.linalg.norm(gf_tau - gf_tau_fit, ord=2) ** 2
    return rho_prime.value, gf_tau_fit, chi2

def _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies):
    energies_extract = np.linspace(emin, emax, num=num_energies)
    delta_energy = energies_extract[1] - energies_extract[0]
    kernel = _get_kernel_matrix(energies_extract, tau_grid, beta, delta_energy)
    U, S, Vt = _getSVD(kernel, nsv=nsv)
    return U, S, Vt, delta_energy, energies_extract

def _integral_kramers_kronig(energies_imagpart, energy_realpart, gf_imag_interp, energy_threshold):
    enum1 = gf_imag_interp(energies_imagpart) - gf_imag_interp(energy_realpart)
    enum2 = np.gradient(gf_imag_interp(energies_imagpart), energies_imagpart)
    den = energies_imagpart - energy_realpart
    mask = np.where(np.abs(den) < energy_threshold, 1, 0) #1 if den is below threshold
    #avoid divide by zero with mask * energy_threshold, in this case (1 - mask == 0)
    term1 = enum1 / (den + mask * energy_threshold) * (1 - mask)
    term2 = enum2 * mask 
    kernel = term1 + term2
    integral = np.trapz(y=kernel, x=energies_imagpart)
    return integral

def get_kramers_kronig_realpart(energies, gf_imag, energy_threshold=1e-10, dos_threshold=1e-5):
    if np.abs(gf_imag[0]) > dos_threshold:
        print('Warning! DOS at left interval end exceeds {}.'.format(dos_threshold))
    if np.abs(gf_imag[-1]) > dos_threshold:
        print('Warning! DOS at right interval end exceeds {}.'.format(dos_threshold))
    gf_real = np.zeros(energies.shape, dtype=np.float64)
    gf_imag_interp = interp1d(energies, gf_imag, kind='linear', bounds_error=False, fill_value=0.0)
    energies_noend = energies[1:-1]
    a, b = energies[0], energies[-1]
    gf_real[1:-1] = gf_imag_interp(energies_noend) / np.pi * np.log((b - energies_noend) / (energies_noend - a)) #assume zero DOS at endpoints
    integral_func = lambda y : _integral_kramers_kronig(energies, y, gf_imag_interp, energy_threshold)
    gf_real += np.vectorize(integral_func)(energies) / np.pi
    return energies, gf_real, gf_imag

def dos_to_gf_imag(dos):
    return -np.pi * dos
