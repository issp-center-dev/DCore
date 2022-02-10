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


import numpy
from scipy.linalg import block_diag
import os
import shutil
from itertools import product

from dcore._dispatcher import *

from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0, get_block_size, float_to_complex_array, make_hermite_conjugate, expand_path
from .base import SolverBase


def remove_positive_eigenvalues(Delta_tau):
    ntau = Delta_tau.shape[0]

    for itau in range(ntau):
        evals, evecs = numpy.linalg.eigh(Delta_tau[itau, :, :])
        evals[evals>0] = 0.0
        Delta_tau[itau, :, :] = evecs.dot(numpy.diag(evals).dot(evecs.transpose().conjugate()))

def to_numpy_array(g, block_names):
    """
    Convert BlockGf object to numpy.
    Rearrange spins and orbitals so that up and down spins appear alternatingly.
    If there is a single block, we assume that spin and down spins appear alternatignly.
    If there are two blocks, we assume that they are spin1 and spin2 sectors.
    
    Note by HS: The comment above may be wrong. The correct one may be
       If there is a single block, we assume that a up-spin block is followed by a down-spin block.
    """

    if g.n_blocks > 2:
        raise RuntimeError("n_blocks must be 1 or 2.")

    block_sizes = [get_block_size(g[name]) for name in block_names]
    n_spin_orbital = numpy.sum(block_sizes)

    # FIXME: Bit ugly
    n_data = g[block_names[0]].data.shape[0]

    data = numpy.zeros((n_data, n_spin_orbital, n_spin_orbital), dtype=complex)
    offset = 0
    for ib, name in enumerate(block_names):
        block = g[name]
        block_dim = block_sizes[ib]
        data[:, offset:offset + block_dim, offset:offset + block_dim] = block.data
        offset += block_dim

    # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
    # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
    norb = int(n_spin_orbital//2)
    index = numpy.zeros((n_spin_orbital), dtype=int)
    index[0::2] = numpy.arange(norb)
    index[1::2] = numpy.arange(norb) + norb
    # Swap cols and rows
    return (data[:, :, index])[:, index, :]


def assign_from_numpy_array_legendre(g, data, block_names):
    """
    Set numpy data to BlockGf of legendre
    """

    if g.n_blocks > 2:
        raise RuntimeError("n_blocks must be 1 or 2.")

    block_sizes = [get_block_size(g[name]) for name in block_names]
    n_spin_orbital = numpy.sum(block_sizes)

    assert data.shape[-1] == g[block_names[0]].data.shape[0]

    norb = n_spin_orbital//2
    index = numpy.zeros((n_spin_orbital), dtype=int)
    index[:norb] = 2 * numpy.arange(norb)
    index[norb:] = 2 * numpy.arange(norb) + 1
    # * Change the order of spins and orbitals
    # * Move the last index for $l$ to the first index.
    data_rearranged = (data[:, index, :][index, :, :]).transpose((2,0,1))

    offset = 0
    for name in block_names:
        block = g[name]
        block_dim = get_block_size(block)
        block.data[:,:,:] = data_rearranged[:, offset:offset + block_dim, offset:offset + block_dim]
        offset += block_dim


class ALPSCTHYBSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.

        """

        super(ALPSCTHYBSolver, self).__init__(beta, gf_struct, u_mat, n_iw)
        self.n_tau = max(10001, 10 * n_iw)

    def solve(self, rot, mpirun_command, params_kw):
        """
        In addition to the parameters described in the docstring of SolverBase,
        params_kw must may contain the following parameters.
          exec_path : str, path to an executable, mandatory
          dry_run   : bool, actual computation is not performed if dry_run is True, optional
        """

        self._solve_impl(rot, mpirun_command, None, params_kw)

    def calc_Xloc_ph(self, rot, mpirun_command, num_wf, num_wb, params_kw):
        raise RuntimeError("calc_Xloc_ph is not implemented!")

    def calc_G2loc_ph_sparse(self, rot, mpirun_command, wsample_ph, params_kw):
        self._solve_impl(rot, mpirun_command, wsample_ph, params_kw)

        return self._G2loc_ph_sparse

    def _solve_impl(self, rot, mpirun_command, wsample_ph, params_kw):
        """
        In addition to the parameters described in the docstring of SolverBase,
        params_kw must may contain the following parameters.
          exec_path : str, path to an executable, mandatory
          dry_run   : bool, actual computation is not performed if dry_run is True, optional
        """

        internal_params = {
            'exec_path'           : '',
            'random_seed_offset'  : 0,
            'dry_run'             : False,
        }

        def _read(key):
            if key in params_kw:
                return params_kw[key]
            else:
                return internal_params[key]

        if 'max_time' in params_kw:
            raise RuntimeError("Parameter max_time has been replaced by timelimit!")

        if not 'timelimit' in params_kw:
            raise RuntimeError("Please set timelimit for ALPS/cthyb!")

        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin2, spin1, spin2, ...
        H0 = extract_H0(self._G0_iw, self.block_names)

        # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
        # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
        index = numpy.zeros((2*self.n_orb), dtype=int)
        index[0::2] = numpy.arange(self.n_orb)
        index[1::2] = numpy.arange(self.n_orb) + self.n_orb
        # Swap cols and rows
        H0 = (H0[:, index])[index, :]

        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n + mu - H0 - G0^{-1}(iw_n)
        # H0 -mu is extracted from the tail of G0. (correct?)
        self._Delta_iw = delta(self._G0_iw)
        Delta_tau = make_block_gf(GfImTime, self.gf_struct, self.beta, self.n_tau)
        for name in self.block_names:
            Delta_tau[name] << Fourier(self._Delta_iw[name])
        Delta_tau_data = to_numpy_array(Delta_tau, self.block_names)
        remove_positive_eigenvalues(Delta_tau_data)

        # non-zero elements of U matrix
        # Note: notation differences between ALPS/CT-HYB and TRIQS!
        #    The positions of l and k are swapped.
        U_nonzeros = []
        def conv(i):
            if i < self.n_orb:
                return 2*i
            else:
                return 2*(i-self.n_orb) + 1
        for i, j, k, l in product(list(range(self.n_flavors)), repeat=4):
            if numpy.abs(self.u_mat[i,j,k,l]) > 1e-10:
                alps_index = (conv(i), conv(j), conv(l), conv(k)) # Here, l and k are swapped.
                U_nonzeros.append((alps_index, self.u_mat[i,j,k,l]))

        if rot is None:
            rot_mat_alps = numpy.identity(2*self.n_orb, dtype=complex)
        else:
            assert isinstance(rot, dict)
            if self.use_spin_orbit:
                assert len(rot) == 1
                rot_single_block = rot['ud']
            else:
                rot_single_block = block_diag(*[rot[name] for name in self.block_names])
            rot_mat_alps = numpy.zeros((2*self.n_orb, 2*self.n_orb), dtype=complex)
            for i, iv in product(list(range(2*self.n_orb)), repeat=2):
                rot_mat_alps[conv(i), iv] = rot_single_block[i, iv]

        # Set up input parameters for ALPS/CT-HYB
        p_run = {
            'SEED'                            : params_kw['random_seed_offset'],
            'model.sites'                     : self.n_orb,
            'model.spins'                     : 2,
            'model.beta'                      : self.beta,
            'model.hopping_matrix_input_file' : './hopping.txt',
            'model.coulomb_tensor_input_file' : './Uijkl.txt',
            'model.basis_input_file'          : './basis.txt',
            'model.n_tau_hyb'                 : self.n_tau - 1,
            'model.delta_input_file'          : './delta.txt',
            'measurement.G1.n_tau'            : self.n_tau - 1,
            'measurement.G1.n_matsubara'      : self.n_iw,
        }

        if not wsample_ph is None:
            p_run['measurement.G2.matsubara.on'] = 1
            p_run['measurement.G2.matsubara.frequencies_PH'] = './freqs_PH.txt'

        if os.path.exists('./input.out.h5'):
            shutil.move('./input.out.h5', './input_prev.out.h5')

        # Set parameters specified by the user
        for k, v in list(params_kw.items()):
            if k in internal_params:
                continue
            if k in p_run:
                raise RuntimeError("Cannot override input parameter for ALPS/CT-HYB: " + k)
            p_run[k] = v

        with open('./input.ini', 'w') as f:
            for k, v in list(p_run.items()):
                print(k, " = ", v, file=f)

        with open('./hopping.txt', 'w') as f:
            for i, j in product(list(range(self.n_flavors)), repeat=2):
                print('{} {} {:.15e} {:.15e}'.format(i, j, H0[i,j].real, H0[i,j].imag), file=f)

        with open('./delta.txt', 'w') as f:
            for itau, f1, f2 in product(list(range(self.n_tau)), list(range(self.n_flavors)), list(range(self.n_flavors))):
                print('{} {} {} {:.15e} {:.15e}'.format(itau, f1, f2, Delta_tau_data[itau, f1, f2].real, Delta_tau_data[itau, f1, f2].imag), file=f)

        with open('./Uijkl.txt', 'w') as f:
            print(len(U_nonzeros), file=f)
            for n, elem in enumerate(U_nonzeros):
                i, j, k, l = elem[0]
                print('{} {} {} {} {} {:.15e} {:.15e}'.format(n, i, j, k, l, elem[1].real, elem[1].imag), file=f)

        with open('./basis.txt', 'w') as f:
            for f1, f2 in product(list(range(self.n_flavors)), list(range(self.n_flavors))):
                print('{} {} {:.15e} {:.15e}'.format(f1, f2, rot_mat_alps[f1, f2].real, rot_mat_alps[f1, f2].imag), file=f)

        if not wsample_ph is None:
            nsmpl_freq = wsample_ph[0].shape[0]
            with open('./freqs_PH.txt', 'w') as f:
                print(nsmpl_freq, file=f)
                for i in range(nsmpl_freq):
                    print('{}  {} {} {}'.format(i, wsample_ph[0][i]//2, wsample_ph[1][i]//2, wsample_ph[2][i]//2), file=f)

        if _read('dry_run'):
            return

        # Invoke subprocess
        exec_path = expand_path(_read('exec_path'))

        # Run a working horse
        with open('./output', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, './input.ini'], output_f)

        with open('./output', 'r') as output_f:
            for line in output_f:
                print(line, end='')

        # Read the computed Green's function in Legendre basis and compute G(iwn)
        if not os.path.exists('./input.out.h5'):
            raise RuntimeError("Output HDF5 file of ALPS/CT-HYB does not exist. Something went wrong!")

        with HDFArchive('input.out.h5', 'r') as f:
            # Sign
            sign = f['Sign']
            print("Average sign is ", sign, ".")
            if numpy.abs(sign) < 0.01:
                print("Average sign may be too small!")

            # G_l in Legendre basis
            gl_data = float_to_complex_array(f['G1_LEGENDRE'])
            Nl = gl_data.shape[2]
            G_l = make_block_gf(GfLegendre, self.gf_struct, self.beta, n_points=Nl)
            assign_from_numpy_array_legendre(G_l, gl_data, self.block_names)

            # G_iw with 1/iwn tail
            for name in self.block_names:
                self._Gimp_iw[name] << LegendreToMatsubara(G_l[name])
            make_hermite_conjugate(self._Gimp_iw)

            # Two-particle GF
            if not wsample_ph is None:
                self._G2loc_ph_sparse = (f['G2']['matsubara']['data']).\
                    reshape((self.n_flavors, self.n_flavors, self.n_flavors, self.n_flavors, -1))

        # Solve Dyson's eq to obtain Sigma_iw
        # Sigma_iw = G0_iw^{-1} - G_imp_iw^{-1}
        self._Sigma_iw = dyson(G0_iw=self._G0_iw, G_iw=self._Gimp_iw)


    def name(self):
        return "ALPS/cthyb"


    def get_Delta_iw(self):
        return self._Delta_iw.copy()
