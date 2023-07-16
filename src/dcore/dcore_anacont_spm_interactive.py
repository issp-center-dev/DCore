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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from dcore._dispatcher import MeshReFreq, MeshImFreq, GfReFreq, GfImFreq
from dcore.version import version, print_header
from dcore.anacont_spm import calc_gf_tau_from_gf_matsubara, get_single_continuation, get_multiple_continuations, get_kramers_kronig_realpart, dos_to_gf_imag

def _plot_overview(lambdas, chi2_values, energies, densities, nrows=3, ncols=5):
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, figsize=(15, 10))
    
    axes[0, 0].plot(np.log10(lambdas), chi2_values, '-o', color='k')
    axes[0, 0].set_yscale('log')
    axes[0, 0].set_xlabel(r'$\log_{10} \lambda$')
    axes[0, 0].set_ylabel(r'$\log_{10} \chi^2 $')
    axes[0, 0].xaxis.set_minor_locator(MultipleLocator(1))

    for i, (d, lambd) in enumerate(zip(densities, np.log10(lambdas)), start=1):
        j = i // ncols
        k = i % ncols
        ax = axes[j, k]
        ax.plot(energies, d)
        ax.set_xlim(energies[0], energies[-1])
        ax.set_ylim(0, None)
        ax.set_xlabel(r'$\omega$')
        ax.set_ylabel(r'$\rho(\omega)$')
        ax.text(0.98, 0.87, r'$\log_{{10}} \lambda = {0:1.1f}$'.format(lambd), horizontalalignment='right', transform = ax.transAxes)

    plt.tight_layout()
    plt.show()
    plt.close()

def _anacont_spm_per_gf(params, matsubara_frequencies, gf_matsubara):
    tau_grid, gf_tau, const_real_tail, const_imag_tail = calc_gf_tau_from_gf_matsubara(matsubara_frequencies, gf_matsubara, params['spm_interactive']['n_tau'], params['spm_interactive']['n_tail'], params['beta'], show_fit=params['spm_interactive']['show_fit'])
    num_lambdas = params['spm_interactive']['n_rows_overview'] * params['spm_interactive']['n_cols_overview'] - 1
    lambda_values = np.logspace(start=params['spm_interactive']['lambda_min_log10'], stop=params['spm_interactive']['lambda_max_log10'], num=num_lambdas, base=10, endpoint=True)
    energies, densities, chi2_values, rho_integrated_values = get_multiple_continuations(tau_grid, gf_tau, params['spm_interactive']['n_sv'], params['beta'], params['omega_min'], params['omega_max'], params['Nomega'], const_imag_tail, lambda_values, verbose=params['spm_interactive']['verbose_opt'], max_iters=params['spm_interactive']['max_iters_opt'], solver=params['spm_interactive']['solver_opt'])

    _plot_overview(lambda_values, chi2_values, energies, densities, nrows=params['spm_interactive']['n_rows_overview'], ncols=params['spm_interactive']['n_cols_overview'])
    lambd = input('\nPlease enter desired value for log10(lambda): ')
    lambd = float(lambd)
    print(f'Calculating result with log10(lambda)={lambd}...')
    lambd = 10 ** lambd

    density, _, energies, _, _ = get_single_continuation(tau_grid, gf_tau, params['spm_interactive']['n_sv'], params['beta'], params['omega_min'], params['omega_max'], params['Nomega'], const_imag_tail, lambd, verbose=params['spm_interactive']['verbose_opt'], max_iters=params['spm_interactive']['max_iters_opt'], solver=params['spm_interactive']['solver_opt'])
    
    energies, gf_real, gf_imag = get_kramers_kronig_realpart(energies, dos_to_gf_imag(density))
    gf_real += const_real_tail
    return energies, gf_real, gf_imag

def dcore_anacont_spm_interactive(seedname):
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
            print(f'Performing analytic continuation for data index {idata} and orbital index {i_orb}...')
            gf_matsubara = gf_iw.data[n_matsubara:, i_orb, i_orb]
            plt.plot(matsubara_frequencies, -np.imag(gf_matsubara))
            plt.title('Please determine energy up to which high-energy tail behaves linearly.')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$ \omega_n $')
            plt.ylabel(r'-Im $G( \omega_n )$')
            plt.show()
            plt.close()
            tail_cutoff_energy = input('\nPlease enter energy up to which high-energy tail behaves linearly: ')
            tail_cutoff_energy = float(tail_cutoff_energy)
            i_cutoff_energy = np.searchsorted(matsubara_frequencies, tail_cutoff_energy)
            print('Retaining {} Matsubara frequencies.'.format(i_cutoff_energy))
            retained_matsubara_frequencies = matsubara_frequencies[:i_cutoff_energy]
            gf_matsubara = gf_matsubara[:i_cutoff_energy]
            energies, gf_real, gf_imag = _anacont_spm_per_gf(params, retained_matsubara_frequencies, gf_matsubara)
            sigma_w_data[i_orb, i_orb, :] = gf_real + 1j * gf_imag
            if params['spm_interactive']['show_result']:
                plt.axhline(y=0, xmin=energies[0], xmax=energies[-1], color='lightgrey')
                plt.plot(energies, gf_real, label=r'Re $G( \omega )$')
                plt.plot(energies, gf_imag, label=r'Im $G( \omega )$')
                plt.xlim(energies[0], energies[-1])
                plt.xlabel(r'$\omega$')
                plt.legend()
                plt.tight_layout()
                plt.show()
                plt.close()
        sigma_w = GfReFreq(mesh=mesh_w, data=sigma_w_data)
        data_w[key] = sigma_w.data

    print('Writing to', seedname + '_sigma_w.npz...')
    np.savez(seedname + '_sigma_w.npz', **data_w)

#example file for 'seedname_anacont.toml'
# beta = 40.0
# Nomega = 4000
# omega_min = -6.0
# omega_max = 6.2

# [spm_interactive]
# n_tail = 100
# n_tau = 10000
# n_sv = 30
# show_fit = false
# show_result = true
# lambda_min_log10 = -10
# lambda_max_log10 = 2
# verbose_opt = false
# max_iters_opt = 100
# solver_opt = 'ECOS'
# n_rows_overview = 3
# n_cols_overview = 5

def run():
    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_anacont_spm_interactive.py',
        description='pre script for dcore.',
        usage='$ dcore_anacont_spm_interactive input',
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

    dcore_anacont_spm_interactive(args.seedname)
