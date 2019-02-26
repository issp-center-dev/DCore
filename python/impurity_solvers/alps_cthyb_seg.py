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
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0
from .base import SolverBase


def to_numpy_array(g):
    """
    Convert BlockGf object to numpy.
    Rearrange spins and orbitals so that up and down spins appear alternatively.
    If there is a single block, we assume that spin and down spins appear alternatively.
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
    index = numpy.zeros(n_spin_orbital, dtype=int)
    index[0::2] = numpy.arange(norb)
    index[1::2] = numpy.arange(norb) + norb
    # Swap cols and rows
    return (data[:, :, index])[:, index, :]


def assign_from_numpy_array(g, data):
    """
    Does inversion of to_numpy_array
    data[spin,orb,iw]
    g[spin].data[iw,orb1,orb2] 
    """
    print(g.n_blocks)
    if g.n_blocks != 2:
        raise RuntimeError("n_blocks must be 1 or 2.")

    names = [name for name, block in g]
    norb = data.shape[1]
    niw = data.shape[2]
    print(data.shape)
    print("norb:", norb)

    #check number of Matsubara frequency
    assert data.shape[2]*2 == g[names[0]].data.shape[0]
    print(g[names[0]].data.shape)

    for spin in range(2):
        for orb in range(norb):
            print(orb, spin, names[spin])
            #positive frequency
            g[names[spin]].data[niw:, orb, orb] = data[spin][orb][:]
            #negative frequency
            g[names[spin]].data[:niw, orb, orb] = numpy.conj(data[spin][orb][::-1])


def dcore2alpscore(dcore_U):

    dcore_U_len = len(dcore_U)
    alps_U  = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_Uprime = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_J = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)

    #m_range = range(size)
    for i, j in product(range(dcore_U_len), range(dcore_U_len)):
        alps_U[i, j] = dcore_U[i, j, i, j].real - dcore_U[i, j, j, i].real
        alps_Uprime[i, j] = dcore_U[i, j, i, j].real
        alps_J[i, j] = dcore_U[i, j, j, i].real

    return alps_U, alps_Uprime, alps_J


def write_Umatrix(U, Uprime, J, norb):
    Uout = numpy.zeros((2, norb, 2, norb))
    for a1, a2 in product(range(norb), repeat=2):
        for s1, s2 in product(range(2), repeat=2): #spin-1/2
            if s1 == s2:
                Uout[s1, a1, s2, a2] = U[2*a1+s1, 2*a2+s2]
            else:
                Uout[s1, a1, s2, a2] = Uprime[2*a1+s1, 2*a2+s2] - J[2*a1+s1, 2*a2+s2]

    Uout = Uout.reshape((2*norb, 2*norb))
    with open('./Umatrix', 'w') as f:
        for i in range(2*norb):
            for j in range(2*norb):
                print('{:.15e} '.format(Uout[i, j].real), file=f, end="")
            print("", file=f)


class ALPSCTHYBSEGSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025, n_tau=10001):
        """
        Initialize the solver.
        """

        super(ALPSCTHYBSEGSolver, self).__init__(beta, gf_struct, u_mat, n_iw, n_tau)

    def solve(self, rot, mpirun_command, params_kw):
        """
        In addition to the parameters described in the docstring of SolverBase,
        one can pass solver-dependent parameters using params_kw. For example,
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
        H0 = extract_H0(self._G0_iw)

        # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
        # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
        index = numpy.zeros((2*self.n_orb), dtype=int)
        index[0::2] = numpy.arange(self.n_orb)
        index[1::2] = numpy.arange(self.n_orb) + self.n_orb
        # Swap cols and rows
        H0 = (H0[:, index])[index, :]

        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)
        # H0 is extracted from the tail of the Green's function.
        self._Delta_iw = delta(self._G0_iw)
        Delta_tau = make_block_gf(GfImTime, self.gf_struct, self.beta, self.n_tau)
        for name, block in self._Delta_iw:
            Delta_tau[name] << InverseFourier(self._Delta_iw[name])
        Delta_tau_data = to_numpy_array(Delta_tau)

        # (1c) Set U_{ijkl} for the solver
        # Set up input parameters and files for ALPS/CTHYB-SEG

        p_run = {
            'SEED'                            : params_kw['random_seed_offset'],
            'FLAVORS'                         : self.n_orb*2,
            'BETA'                            : self.beta,
            'N'                               : self.n_tau - 1,
            'NMATSUBARA'                      : self.n_iw,
            'U_MATRIX'                        : 'Umatrix',
            'MU_VECTOR'                       : 'MUvector',
            'cthyb.DELTA'                     : 'delta',
        }

        if os.path.exists('./input.out.h5'):
            shutil.move('./input.out.h5', './input_prev.out.h5')
        # Set parameters specified by the user
        for k, v in params_kw.items():
            if k in internal_params:
                continue
            if k in p_run:
                raise RuntimeError("Cannot override input parameter for ALPS/CTHYB-SEG: " + k)
            p_run[k] = v

        with open('./input.ini', 'w') as f:
            for k, v in p_run.items():
                print(k, " = ", v, file=f)

        with open('./delta', 'w') as f:
            for itau in range(self.n_tau):
                print('{}'.format(itau), file=f, end="")
                for f1 in range(self.n_flavors):
                    if Delta_tau_data[itau, f1, f1].real >0:
                        Delta_tau_data[itau, f1, f1] = 0
                    print(' {:.15e}'.format(Delta_tau_data[itau, f1, f1].real), file=f, end="")
                print("", file=f)

        U, Uprime, J = dcore2alpscore(self.u_mat)
        write_Umatrix(U, Uprime, J, self.n_orb)
        
        with open('./MUvector', 'w') as f:
            for orb in range(self.n_orb):
                for spin in range(2):
                    print('{:.15e} '.format(H0[2*orb+spin][2*orb+spin].real), file=f, end="")
            print("", file=f)

        if _read('dry_run'):
            return

        # Invoke subprocess
        exec_path = _read('exec_path')
        if exec_path == '':
            raise RuntimeError("Please set exec_path!")
        if not os.path.exists(exec_path):
            raise RuntimeError(exec_path + " does not exist. Set exec_path properly!")

        # (2) Run a working horse
        with open('./output', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, './input.ini'], output_f)

        with open('./output', 'r') as output_f:
            for line in output_f:
                print(line, end='')

        # (3) Copy results into
        #   self._Sigma_iw
        swdata = numpy.zeros((2, self.n_orb, self.n_iw), dtype=complex)
        with HDFArchive('sim.h5', 'r') as f:
            for orb in range(self.n_orb):
                for spin in range(2):
                    swdata_array = f["S_omega"][str(orb*2+spin)]["mean"]["value"]
                    assert swdata_array.dtype == numpy.complex
                    assert swdata_array.shape == (self.n_iw,)
                    swdata[spin, orb, :] = swdata_array
        assign_from_numpy_array(self._Sigma_iw, swdata)

    def name(self):
        return "ALPS/cthyb-seg"
