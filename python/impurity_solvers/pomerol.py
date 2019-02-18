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
from itertools import product
import os

from pytriqs.gf.local import *
# from pytriqs.archive import HDFArchive
from pytriqs.operators import *

from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0
from .base import SolverBase

# bse
from dft_tools.index_pair import IndexPair, IndexPair2
from bse_tools.h5bse import h5BSE

def assign_from_numpy_array(g_block, data):

    # g_block[].data.shape = (2*n_iw, n_orb, n_orb)
    #                        (2*n_iw, 2*n_orb, 2*n_orb)
    # Both positive and negative freq are stored

    # data.shape = (2, n_orb, n_orb, n_iw)       w/ spin-orbit
    #              (1, 2*n_orb, 2*n_orb, n_iw)   w/o spin-orbit
    # Only positive freq

    for i, (bname, gf) in enumerate(g_block):
        # FIXME: spin order

        # print(bname)
        # print(gf.data.shape)
        # print(data[i].shape)

        # number of positive Matsubara freq
        n_iw = data[i].shape[2]
        assert gf.data.shape[0] == 2*n_iw

        # data(o1, o2, iw) --> (iw, o1, o2)
        gf_iw_o1_o2 = data[i].transpose(2, 0, 1)
        # print(gf_iw_o1_o2.shape)

        # copy data in the positive freq side
        gf.data[n_iw:, :, :] = gf_iw_o1_o2.copy()

        # copy data in the negative freq side (complex conjugate)
        gf.data[0:n_iw, :, :] = gf_iw_o1_o2.transpose(0, 2, 1)[::-1, :, :].conjugate().copy()


def set_tail(g_block):
    # TODO: compute tails in ED
    for bname, gf in g_block:
        gf.tail.zero()
        gf.tail[1] = numpy.identity(gf.N1)


def decompose_index(index, n_orb):
    spn = index / n_orb
    orb = index % n_orb
    print(spn, orb)
    return spn, orb

class PomerolSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025, n_tau=10001):
        """
        Initialize the solver.
        """

        super(PomerolSolver, self).__init__(beta, gf_struct, u_mat, n_iw, n_tau)

    def solve(self, rot, mpirun_command, params_kw):
        """
        In addition to the parameters described in the docstring of SolverBase,
        one can pass solver-dependent parameters using params_kw. For example,
          exec_path : str, path to an executable, mandatory
          dry_run   : bool, actual computation is not performed if dry_run is True, optional
        """

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

        # print("params_kw =", params_kw)
        exec_path = params_kw['exec_path']
        flag_vx = params_kw.get('flag_vx', 0)
        dir_vx = params_kw.get('dir_vx', './two_particle')
        n_w2f = params_kw.get('n_w2f', 10)
        n_w2b = params_kw.get('n_w2b', 1)

        # for BSE
        # n_corr_shells = params_kw.get('n_corr_shells', None)
        # icrsh = params_kw.get('icrsh', None)
        # only_diagonal = not params_kw.get('nonlocal_order_parameter', False)
        # FIXME: passed from DCore
        n_corr_shells = params_kw.get('n_corr_shells', 1)
        icrsh = params_kw.get('icrsh', 0)
        only_diagonal = not params_kw.get('nonlocal_order_parameter', False)

        # print(self.gf_struct)
        flag_spin_conserve = 1 if len(self.gf_struct) == 2 else 0

        file_pomerol = "pomerol.in"
        file_h0 = "h0.in"
        file_umat = "umat.in"
        file_gf = "gf.dat"

        params_pomerol = {}
        params_pomerol['n_orb'] = self.n_orb
        params_pomerol['beta'] = self.beta
        params_pomerol['flag_spin_conserve'] = flag_spin_conserve
        params_pomerol['file_h0'] = file_h0
        params_pomerol['file_umat'] = file_umat
        params_pomerol['flag_gf'] = 1
        params_pomerol['n_w'] = self.n_iw
        params_pomerol['flag_vx'] = flag_vx
        params_pomerol['n_w2f'] = n_w2f
        params_pomerol['n_w2b'] = n_w2b

        with open(file_pomerol, "w") as f:
            for key, val in params_pomerol.items():
                print(key, "=", val, file=f)

        # (1a) If H0 is necessary:
        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin1, ..., spin2, spin2, ...
        H0 = extract_H0(self._G0_iw)
        H0 = 0.5 * (H0.transpose().conj() + H0)
        # print(H0.shape)
        # print(H0)
        with open(file_h0, "w") as f:
            for i, j in product(range(H0.shape[0]), range(H0.shape[1])):
                # TODO: check order of spin/orbital
                # TODO: real or complex
                if abs(H0[i,j]) != 0:
                    # print(i, j, H0[i,j].real, H0[i,j].imag, file=f)
                    print(i, j, H0[i,j].real, file=f)

        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)
        # H0 is extracted from the tail of the Green's function.
        self._Delta_iw = delta(self._G0_iw)

        # (1c) Set U_{ijkl} for the solver

        # print(self.u_mat)
        with open(file_umat, "w") as f:
            for i, j, k, l in product(range(self.n_flavors), repeat=4):
                # TODO: check order of spin/orbital
                # TODO: real or complex
                if abs(self.u_mat[i,j,k,l]) != 0:
                    # print(i, j, k, l, self.u_mat[i, j, k, l].real, self.u_mat[i, j, k, l].imag, file=f)
                    print(i, j, k, l, self.u_mat[i, j, k, l].real, file=f)

        # (2) Run a working horse
        with open('./stdout.log', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, file_pomerol], output_f)

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

        # load data as a complex type
        gf_1d = numpy.loadtxt(file_gf).view(complex).reshape(-1)
        # print(gf_1d.shape)
        if flag_spin_conserve:
            gf = gf_1d.reshape((2, self.n_orb, self.n_orb, self.n_iw))
        else:
            gf = gf_1d.reshape((1, self.n_flavors, self.n_flavors, self.n_iw))
        # print(gf.shape)
        assign_from_numpy_array(self._Gimp_iw, gf)

        set_tail(self._Gimp_iw)

        # Compute Sigma_iw
        # NOTE: do NOT use self._G0_iw because it includes more information than that passed to the solver
        n_block = len(self.gf_struct)
        n_inner = H0.shape[0] / n_block
        # cut H0 into block structure
        h0_block = [ H0[s*n_inner:(s+1)*n_inner, s*n_inner:(s+1)*n_inner] for s in range(n_block) ]
        g0_inv = make_block_gf(GfImFreq, self.gf_struct, self.beta, self.n_iw)
        g0_inv << iOmega_n
        g0_inv -= h0_block
        self._Sigma_iw << g0_inv - inverse(self._Gimp_iw)

        # Solve Dyson's eq to obtain Sigma_iw
        # self._Sigma_iw = dyson(G0_iw=self._G0_iw, G_iw=self._Gimp_iw)

        if flag_vx:
            block2 = IndexPair2(range(n_corr_shells), range(2), only_diagonal1=only_diagonal)
            inner2 = IndexPair(range(self.n_orb), convert_to_int=True)

            xloc = numpy.zeros((self.n_orb**2, self.n_orb**2, 2*n_w2f, 2*n_w2f))

            # TODO: spin loop first
            for i1,i2,i3,i4 in product(range(self.n_flavors), repeat=4):
                filename = dir_vx + "/%d_%d_%d_%d.dat" %(i1,i2,i3,i4)
                # print(filename)
                if os.path.exists(filename):
                    print(filename)
                    # load data as a complex type
                    data = numpy.loadtxt(filename).view(complex).reshape(-1)
                    print(data.shape)
                    data = data.reshape((n_w2b, 2*n_w2f, 2*n_w2f))
                    print(data.shape)

                    s1, o1 = decompose_index(i1, self.n_orb)
                    s2, o2 = decompose_index(i2, self.n_orb)
                    s3, o3 = decompose_index(i3, self.n_orb)
                    s4, o4 = decompose_index(i4, self.n_orb)

                    s12 = block2.get_index(icrsh, s1, icrsh, s2)
                    s43 = block2.get_index(icrsh, s4, icrsh, s3)
                    inner12 = inner2.get_index(o1, o2)
                    inner43 = inner2.get_index(o4, o3)

                    wb = 0
                    xloc[inner12, inner43, :, :] = data[wb, :, :]
                    # X_loc[(s12, s43)] =

                    # h5bse
                    # BS = h5BSE(h5_file, groupname)
                    # BS.save(key=('X_loc', wb), data=X_loc)


    def name(self):
        return "pomerol"

    def get_Delta_iw(self):
        return self._Delta_iw.copy()
