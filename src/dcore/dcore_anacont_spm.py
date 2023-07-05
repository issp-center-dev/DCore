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

import argparse
import toml
from dcore._dispatcher import MeshReFreq, MeshImFreq, GfReFreq, GfImFreq
from dcore.version import version, print_header
import numpy as np
import cvxpy as cp
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from scipy.sparse.linalg import svds

def _find_sum_rule_const(matsubara_frequencies, gf_wn, ntail, c0=1.0, show_fit=False):
    wntail = matsubara_frequencies[-ntail:]
    x = 1.0 / wntail
    y = np.imag(gf_wn[-ntail:])
    result = least_squares(lambda c: y + c * x, c0)
    if not result.success:
        print('Finding sum rule constant failed.')
    result = result.x[0]
    if show_fit:
        import matplotlib.pyplot as plt
        z = np.linspace(wntail[0], wntail[-1], num=1000)
        plt.scatter(wntail, y, zorder=5, color='C0')
        plt.plot(z, -result / z, zorder=10, color='C1')
        plt.xlabel(r'$\omega_n$')
        plt.ylabel(r'Im $G( \omega_n )$')
        plt.tight_layout()
        plt.show()
    return result

def _calc_gf_tau(matsubara_frequencies, gf_wn, beta, sum_rule_const, ntau):
    tau_grid = np.linspace(0, beta, num=ntau)
    gf_tau = -0.5 * sum_rule_const * np.ones(tau_grid.shape, dtype=np.float64)
    kernel = lambda tau : np.sum(np.real(gf_wn) * np.cos(matsubara_frequencies * tau) + (np.imag(gf_wn) + sum_rule_const / matsubara_frequencies) * np.sin(matsubara_frequencies * tau))
    gf_tau += 2 / beta * np.vectorize(kernel)(tau_grid)
    gf_tau *= -1
    return tau_grid, gf_tau

def _calc_gf_tau_from_gf_matsubara(matsubara_frequencies, gf_wn, ntau, ntail, beta, show_fit=False):
    sum_rule_const = _find_sum_rule_const(matsubara_frequencies, gf_wn, ntail, show_fit=show_fit)
    print('Determined sum rule constant: {}'.format(sum_rule_const))
    tau_grid, gf_tau = _calc_gf_tau(matsubara_frequencies, gf_wn, beta, sum_rule_const, ntau)
    return tau_grid, gf_tau, sum_rule_const

def _get_single_continuation(tau_grid, gf_tau, nsv, beta, emin, emax, num_energies, sum_rule_const, lambd, verbose=True, max_iters=100, solver='ECOS'):
    U, S, Vt, delta_energy, energies_extract = _get_svd_for_continuation(tau_grid, nsv, beta, emin, emax, num_energies)
    rho_prime, gf_tau_fit, chi2 = _solveProblem(delta_energy, U, S, Vt, gf_tau, sum_rule_const, lambd, verbose=verbose, max_iters=max_iters, solver=solver)
    rho = np.dot(Vt.T, rho_prime)
    rho_integrated = np.trapz(y=rho, x=energies_extract)
    return rho, gf_tau_fit, energies_extract, rho_integrated, chi2

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
    constraints = [Vt.T @ rho_prime >= 0, cp.sum(delta_energy * V_mod @ rho_prime) == sum_rule_const]
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

def _get_kramers_kronig_realpart(energies, gf_imag, energy_threshold=1e-10, dos_threshold=1e-5):
    if np.abs(gf_imag[0]) > dos_threshold:
        print('Warning! DOS at left interval end exceeds {}.'.format(dos_threshold))
    if np.abs(gf_imag[-1]) > dos_threshold:
        print('Warning! DOS at right interval end exceeds {}.'.format(dos_threshold))
    gf_real = np.zeros(energies.shape, dtype=np.float64)
    gf_imag_interp = interp1d(energies, gf_imag, kind='linear', bounds_error=False, fill_value=0.0)
    energies_noend = energies[1:-1]
    a, b = energies[0], energies[-1]
    gf_real[1:-1] = gf_imag_interp(energies_noend) / np.pi * np.log((b - energies_noend) / (energies_noend - a)) #assume zero DOS at endpoints
    integral_func = lambda y : _integral_kramers_kronig(energies, y, gf_imag_interp, energy_threshold) #intentionally use energies for integration grid
    gf_real += np.vectorize(integral_func)(energies) / np.pi
    gf_imag_resampled = gf_imag_interp(energies)
    return energies, gf_real, gf_imag_resampled

def _dos_to_gf_imag(dos):
    return -np.pi * dos

def _anacont_spm_per_gf(params, matsubara_frequencies, gf_matsubara):
    tau_grid, gf_tau, sum_rule_const = _calc_gf_tau_from_gf_matsubara(matsubara_frequencies, gf_matsubara, params['spm']['n_tau'], params['spm']['n_tail'], params['beta'], show_fit=params['spm']['show_fit'])
    density, gf_tau_fit, energies_extract, density_integrated, chi2 = _get_single_continuation(tau_grid, gf_tau, params['spm']['n_sv'], params['beta'], params['omega_min'], params['omega_max'], params['Nomega'], sum_rule_const, params['spm']['lambda'], verbose=params['spm']['verbose_opt'], max_iters=params['spm']['max_iters_opt'], solver=params['spm']['solver_opt'])
    energies, gf_real, gf_imag = _get_kramers_kronig_realpart(energies_extract, _dos_to_gf_imag(density))
    return energies, gf_real, gf_imag

def dcore_anacont_spm(seedname):
    print('Reading ', seedname + '_anacont.toml...')
    with open(seedname + '_anacont.toml', 'r') as f:
        params = toml.load(f)
 
    print('Reading ', seedname + '_sigma_iw.npz...')
    npz = np.load(seedname + '_sigma_iw.npz')

    assert params['omega_min'] < params['omega_max']
    mesh_w = MeshReFreq(params['omega_min'], params['omega_max'], params['Nomega'])

    data_w = {}
    num_data = np.sum([key.startswith('data') for key in npz.keys()])
    for idata in range(num_data):
        key = f'data{idata}'
        data = npz[key]
        n_matsubara = data.shape[0]//2
        mesh_iw = MeshImFreq(params['beta'], 'Fermion', n_matsubara)
        gf_iw = GfImFreq(data=data, beta=params['beta'], mesh=mesh_iw)
        n_orbitals = gf_iw.data.shape[1]
        matsubara_frequencies = np.imag(gf_iw.mesh.values()[n_matsubara:])
        sigma_w_data = np.zeros((n_orbitals, n_orbitals, params['Nomega']), dtype=np.complex128)
        for i_orb in range(n_orbitals):
            print(f'Working on data index {idata} and orbital index {i_orb}...')
            gf_imag_matsubara = gf_iw.data[n_matsubara:, i_orb, i_orb]
            energies, gf_real, gf_imag = _anacont_spm_per_gf(params, matsubara_frequencies, gf_imag_matsubara)
            sigma_w_data[i_orb, i_orb, :] = gf_real + 1j * gf_imag
        sigma_w = GfReFreq(mesh=mesh_w, data=sigma_w_data)
        data_w[key] = sigma_w.data

    print('Writing to', seedname + '_sigma_w.npz...')
    np.savez(seedname + '_sigma_w.npz', **data_w)

def run():
    # install test code with pip3 install .
    # test with DCORE_TRIQS_COMPAT=1 dcore_anacont_spm crgete3
    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_anacont_spm.py',
        description='pre script for dcore.',
        usage='$ dcore_anacont_spm input',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        #epilog=generate_all_description()
    )
    parser.add_argument('seedname',
                        action='store',
                        default=None,
                        type=str,
                        help='seedname'
                        )
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))

    args = parser.parse_args()

    dcore_anacont_spm(args.seedname)
