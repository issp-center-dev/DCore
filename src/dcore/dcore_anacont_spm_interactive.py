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
from dcore._dispatcher import MeshReFreq, MeshImFreq, GfReFreq, GfImFreq
from dcore.version import version, print_header
from dcore.anacont_spm import calc_gf_tau_from_gf_matsubara, get_single_continuation, get_kramers_kronig_realpart, dos_to_gf_imag

def _anacont_spm_per_gf(params, matsubara_frequencies, gf_matsubara):
    tau_grid, gf_tau, sum_rule_const = calc_gf_tau_from_gf_matsubara(matsubara_frequencies, gf_matsubara, params['spm']['n_tau'], params['spm']['n_tail'], params['beta'], show_fit=params['spm']['show_fit'])
    density, gf_tau_fit, energies_extract, density_integrated, chi2 = get_single_continuation(tau_grid, gf_tau, params['spm']['n_sv'], params['beta'], params['omega_min'], params['omega_max'], params['Nomega'], sum_rule_const, params['spm']['lambda'], verbose=params['spm']['verbose_opt'], max_iters=params['spm']['max_iters_opt'], solver=params['spm']['solver_opt'])
    energies, gf_real, gf_imag = get_kramers_kronig_realpart(energies_extract, dos_to_gf_imag(density))
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
            print(f'Working on data index {idata} and orbital index {i_orb}...')
            gf_imag_matsubara = gf_iw.data[n_matsubara:, i_orb, i_orb]
            energies, gf_real, gf_imag = _anacont_spm_per_gf(params, matsubara_frequencies, gf_imag_matsubara)
            sigma_w_data[i_orb, i_orb, :] = gf_real + 1j * gf_imag
            if params['spm']['show_result']:
                import matplotlib.pyplot as plt
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

# [spm]
# n_tail = 100
# n_tau = 10000
# n_sv = 30 # number of retained singular values
# show_fit = false
# show_result = false
# lambda = 1e-6
# verbose_opt = true
# max_iters_opt = 100
# solver_opt = 'ECOS'

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
