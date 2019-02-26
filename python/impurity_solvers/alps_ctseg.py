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
from __future__ import print_function

import numpy
from scipy.linalg import block_diag
import os
from itertools import product

from pytriqs.gf.local import *
# from pytriqs.archive import HDFArchive
from pytriqs.operators import *

from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0
from .base import SolverBase

def to_numpy_array(g):
    """
    Convert BlockGf object to numpy.
    Rearrange spins and orbitals so that up and down spins appear alternatingly.
    If there is a single block, we assume that spin and down spins appear alternatignly.
    If there are two blocks, we assume that they are spin1 and spin2 sectors.
    """

    if g.n_blocks > 2:
        raise RuntimeError("n_blocks must be 1 or 2.")

    names = [name for name, block in g]
    n_spin_orbital = numpy.sum([len(block.indices) for name, block in g])

    # FIXME: Bit ugly
    n_data = g[names[0]].data.shape[0]

    data = numpy.zeros((n_data, n_spin_orbital, n_spin_orbital), dtype=complex)
    offset = 0
    for name, block in g:
        block_dim = len(block.indices)
        data[:, offset:offset + block_dim, offset:offset + block_dim] = block.data
        offset += block_dim

    # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
    # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
    norb = int(n_spin_orbital/2)
    index = numpy.zeros((n_spin_orbital), dtype=int)
    index[0::2] = numpy.arange(norb)
    index[1::2] = numpy.arange(norb) + norb
    # Swap cols and rows
    return (data[:, :, index])[:, index, :]


class ALPSCTSEGSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025, n_tau=10001):
        """
        Initialize the solver.
        """

        super(ALPSCTSEGSolver, self).__init__(beta, gf_struct, u_mat, n_iw, n_tau)

    def solve(self, rot, mpirun_command, params_kw):
        """
        In addition to the parameters described in the docstring of SolverBase,
        one can pass solver-dependent parameters using params_kw. For example,
          exec_path : str, path to an executable, mandatory
          dry_run   : bool, actual computation is not performed if dry_run is True, optional
        """
        internal_params = {
            'exec_path'           : '',
            'work_dir'            : '',
            'random_seed_offset'  : 0,
            'dry_run'             : False,
        }
        def _read(key):
            if key in params_kw:
                return params_kw[key]
            else:
                return internal_params[key]
        print (params_kw)
        # (1) Set configuration for the impurity solver
        # input:
        #   self.beta
        #   self.set_G0_iw
        #   self.u_mat
        #
        # Additionally, the following variables may used:
        #   self.n_orb
        #   self.n_flavor
        #   self.gf_struct
        #   self.n_tau

        # (1a) If H0 is necessary:
        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin1, ..., spin2, spin2, ...
		#TODO: check order spin/orbtital
        H0 = extract_H0(self._G0_iw)
	
        # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
        # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
        index = numpy.zeros((2*self.n_orb), dtype=int)
        index[0::2] = numpy.arange(self.n_orb)
        index[1::2] = numpy.arange(self.n_orb) + self.n_orb
        # Swap cols and rows
        H0 = (H0[:, index])[index, :]

        H0 = 0.5 * (H0.transpose().conj() + H0)

        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
	###### FIX mu?
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)
        # H0 is extracted from the tail of the Green's function.
		#TODO: check order spin/orbtital
        self._Delta_iw = delta(self._G0_iw)
        Delta_tau = make_block_gf(GfImTime, self.gf_struct, self.beta, self.n_tau)
        for name, block in self._Delta_iw:
            Delta_tau[name] << InverseFourier(self._Delta_iw[name])
        Delta_tau_data = to_numpy_array(Delta_tau)

        # non-zero elements of U matrix
        # Note: notation differences between ALPS/CT-HYB and TRIQS!
        #    The positions of l and k are swapped.
        U_nonzeros = []
        def conv(i):
            if i < self.n_orb:
                return 2*i
            else:
                return 2*(i-self.n_orb) + 1
        for i, j, k, l in product(range(self.n_flavors), repeat=4):
            if numpy.abs(self.u_mat[i,j,k,l]) > 1e-10:
                alps_index = (conv(i), conv(j), conv(l), conv(k)) # Here, l and k are swapped.
                U_nonzeros.append((alps_index, self.u_mat[i,j,k,l]))

        if rot is None:
            rot_mat_alps = numpy.identity(2*self.n_orb, dtype=complex)
        else:
            if self.use_spin_orbit:
                rot_single_block = rot
            else:
                rot_single_block = block_diag(*[rot[name] for name in self.block_names])
            rot_mat_alps = numpy.zeros((2*self.n_orb, 2*self.n_orb), dtype=complex)
            for i, j in product(range(2*self.n_orb), repeat=2):
                rot_mat_alps[conv(i), conv(j)] = rot_single_block[i,j]


        # (1c) Set U_{ijkl} for the solver
        # for i, j, k, l in product(range(self.n_flavors), repeat=4):
        #     self.u_mat[i, j, k, l]
        # Set up input parameters for ALPS/CT-HYB
		#TODO: check order spin/orbtital
        #file_umat = "Umatrix.dat"
	#umat_new = numpy.einsum("ijji->ij", self.u_mat)
        #with open(file_umat, "w") as f:
        #    for i, j, k, l in product(range(self.n_flavors), repeat=4):
        #        if abs(umat_new[i,j,k,l]) != 0:
        #            print(umat_new[i, j, k, l].real, file=f)
        
	work_dir = os.path.abspath(params_kw['work_dir'])
        p_run = {
        #MU, Umatrix
            'SEED'                            : params_kw['random_seed_offset'],
            'FLAVORS'                         : self.n_orb*2,
            'BETA'                            : self.beta,
        #debug
         #   'U'                               : params_kw['U'],
         #   'U\''						      : params_kw['U\''],
         #   'cthyb.J'                         : params_kw['J'],
        #debug end
			#FIX currently only works for model calculation
            #'model.hopping_matrix_input_file' : work_dir + '/hopping.txt',
            'N'                               : self.n_tau - 1,
            'NMATSUBARA'                      : self.n_iw,
            't'                               : 1,
            'tprime'                          : 0,
            #fixme
            #'U_MATRIX'                        : work_dir + '/Umatrix.dat',
            #?
            #'model.basis_input_file'          : work_dir + '/basis.txt',
            'cthyb.DELTA'                      : work_dir + '/delta.txt',
        }

        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

            if os.path.exists(work_dir + '/input.out.h5'):
                shutil.move(work_dir + '/input.out.h5', work_dir + '/input_prev.out.h5')
        # Set parameters specified by the user
        for k, v in params_kw.items():
            if k in internal_params:
                continue
            if k in p_run:
                raise RuntimeError("Cannot override input parameter for ALPS/CT-HYB: " + k)
            p_run[k] = v

        with open(work_dir + '/input.ini', 'w') as f:
            for k, v in p_run.items():
                print(k, " = ", v, file=f)

        print("2")
        return
        #with open(work_dir + '/hopping.txt', 'w') as f:
        #    for i, j in product(range(self.n_flavors), repeat=2):
        #        print('{} {} {:.15e} {:.15e}'.format(i, j, H0[i,j].real, H0[i,j].imag), file=f)

        #with open(work_dir + '/delta.txt', 'w') as f:
        #    for itau in product(range(self.n_tau)):
            #for itau, f1, f2 in product(range(self.n_tau), range(self.n_flavors), range(self.n_flavors)):
         #      print('{} {:.15e} {:.15e}'.format(itau, [for f1 in range(self.n_flavors) Delta_tau_data[itau, f1, f1].real]), file=f)
               #print('{} {} {} {:.15e} {:.15e}'.format(itau, f1, f2, Delta_tau_data[itau, f1, f2].real, Delta_tau_data[itau, f1, f2].imag), file=f)

        with open(work_dir + '/Uijkl.txt', 'w') as f:
            print(len(U_nonzeros), file=f)
            for n, elem in enumerate(U_nonzeros):
                i, j, k, l = elem[0]
                print('{} {} {} {} {} {:.15e} {:.15e}'.format(n, i, j, k, l, elem[1].real, elem[1].imag), file=f)

        #with open(work_dir + '/basis.txt', 'w') as f:
        #    for f1, f2 in product(range(self.n_flavors), range(self.n_flavors)):
        #        print('{} {} {:.15e} {:.15e}'.format(f1, f2, rot_mat_alps[f1, f2].real, rot_mat_alps[f1, f2].imag), file=f)

        if _read('dry_run'):
            return

        # Invoke subprocess
        exec_path = _read('exec_path')
        if exec_path == '':
            raise RuntimeError("Please set exec_path!")
        if not os.path.exists(exec_path):
            raise RuntimeError(exec_path + " does not exist. Set exec_path properly!")

        # (2) Run a working horse
        with open(work_dir + '/output', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, work_dir + '/input.ini'], output_f)

        with open(work_dir + '/output', 'r') as output_f:
            for line in output_f:
                print(line, end='')

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

    def name(self):
        return "ALPS/ctseg"

    def get_Delta_iw(self):
        return self._Delta_iw.copy()
