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


def assign_from_numpy_array(g_block, data):

    # DCore:
    #   g_block[].data.shape = (2*n_iw, n_orb, n_orb)
    #                          (2*n_iw, 2*n_orb, 2*n_orb)
    #       Both positive and negative freq are stored

    # pomerol:
    #   data.shape = (2, n_orb, n_orb, n_iw)       w/ spin-orbit
    #                (1, 2*n_orb, 2*n_orb, n_iw)   w/o spin-orbit
    #       Only positive freq

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
        # Additionally, the following variables may be used:
        #   self.n_orb
        #   self.n_flavor
        #   self.gf_struct
        #   self.n_tau
        #   self.use_spin_orbit

        # print("params_kw =", params_kw)
        exec_path = params_kw['exec_path']

        # for BSE
        flag_vx = params_kw.get('flag_vx', 0)
        n_w2f = params_kw.get('num_wf', None)
        n_w2b = params_kw.get('num_wb', None)

        file_pomerol = "pomerol.in"
        file_h0 = "h0.in"
        file_umat = "umat.in"
        file_gf = "gf.dat"

        params_pomerol = {
            'n_orb': self.n_orb,
            'beta': self.beta,
            'flag_spin_conserve': 1 if not self.use_spin_orbit else 0,
            'file_h0': file_h0,
            'file_umat': file_umat,
            'flag_gf': 1,
            'n_w': self.n_iw,
            'flag_vx': flag_vx,
            'n_w2f': n_w2f,
            'n_w2b': n_w2b,
        }

        with open(file_pomerol, "w") as f:
            for key, val in params_pomerol.items():
                print(key, "=", val, file=f)

        # (1a) If H0 is necessary:
        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin1, ..., spin2, spin2, ...
        h0_mat = extract_H0(self._G0_iw)

        with open(file_h0, "w") as f:
            for i, j in product(range(h0_mat.shape[0]), range(h0_mat.shape[1])):
                # TODO: real or complex
                if abs(h0_mat[i, j]) != 0:
                    # print(i, j, H0[i,j].real, H0[i,j].imag, file=f)
                    print(i, j, h0_mat[i, j].real, file=f)

        # (1c) Set U_{ijkl} for the solver
        with open(file_umat, "w") as f:
            for i, j, k, l in product(range(self.n_flavors), repeat=4):
                # TODO: real or complex
                if abs(self.u_mat[i, j, k, l]) != 0:
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
        if not self.use_spin_orbit:
            gf = gf_1d.reshape((2, self.n_orb, self.n_orb, self.n_iw))
        else:
            gf = gf_1d.reshape((1, self.n_flavors, self.n_flavors, self.n_iw))
        assign_from_numpy_array(self._Gimp_iw, gf)

        set_tail(self._Gimp_iw)

        # Compute Sigma_iw
        # NOTE:
        #   compute G0(iw) from h0_mat instead of using self._G0_iw, since
        #   self._G0_iw includes more information than that passed to the
        #   solver (see extract_H0 for details).
        n_block = len(self.gf_struct)
        n_inner = h0_mat.shape[0] / n_block
        # cut H0 into block structure
        h0_block = [h0_mat[s*n_inner:(s+1)*n_inner, s*n_inner:(s+1)*n_inner] for s in range(n_block)]
        g0_inv = make_block_gf(GfImFreq, self.gf_struct, self.beta, self.n_iw)
        g0_inv << iOmega_n
        g0_inv -= h0_block
        self._Sigma_iw << g0_inv - inverse(self._Gimp_iw)

    def calc_g2(self, rot, mpirun_command, params_kw):

        params_kw['flag_vx'] = 1
        assert 'num_wf' in params_kw
        assert 'num_wb' in params_kw

        self.solve(rot, mpirun_command, params_kw)

        from dft_tools.index_pair import IndexPair, IndexPair2
        from bse_tools.h5bse import h5BSE

        n_w2f = params_kw.get('num_wf')
        n_w2b = params_kw.get('num_wb')

        dir_g2 = params_kw.get('dir_g2', './two_particle')
        bse_h5 = params_kw.get('bse_h5_file', 'dmft_bse.h5')
        bse_grp = params_kw.get('bse_grp_name', '')
        bse_info = params_kw.get('bse_info', 'save')

        # FIXME: passed from DCore
        n_corr_shells = params_kw.get('n_corr_shells', 1)
        icrsh = params_kw.get('icrsh', 0)
        only_diagonal = not params_kw.get('nonlocal_order_parameter', False)

        n_spin = 2 if not self.use_spin_orbit else 1
        n_inner = self.n_flavors / n_spin

        # FIXME: spin order and names
        block2 = IndexPair2(range(n_corr_shells), range(n_spin), only_diagonal1=only_diagonal)
        inner2 = IndexPair(range(n_inner), convert_to_int=True)
        print(" block2 namelist =", block2.namelist)
        print(" inner2 namelist =", inner2.namelist)

        # h5bse
        h5bse = h5BSE(bse_h5, bse_grp)
        if bse_info == 'check':
            # check equivalence of info
            assert h5bse.get(key=('block_name',)) == block2.namelist
            assert h5bse.get(key=('inner_name',)) == inner2.namelist
            assert h5bse.get(key=('beta',)) == self.beta
        elif bse_info == 'save':
            # save info
            h5bse.save(key=('block_name',), data=block2.namelist)
            h5bse.save(key=('inner_name',), data=inner2.namelist)
            h5bse.save(key=('beta',), data=self.beta)
        else:
            raise ValueError("bse_info =", bse_info)

        # read X_loc data and save into h5 file
        for wb in range(n_w2b):
            # boson freq
            print(" ---\n wb = %d" % wb)
            x_loc = {}
            for i1, i2, i3, i4 in product(range(self.n_flavors), repeat=4):
                filename = dir_g2 + "/%d_%d_%d_%d.dat" % (i1, i2, i3, i4)
                if not os.path.exists(filename):
                    continue
                print(" reading", filename)

                # load data as a complex type
                # 1d array --> (wb, wf1, wf2) --> (wf1, wf2)
                data = numpy.loadtxt(filename).view(complex).reshape(-1)
                data = data.reshape((n_w2b, 2*n_w2f, 2*n_w2f))
                data_wb = data[wb]

                if not self.use_spin_orbit:
                    s1, o1 = decompose_index(i1, self.n_orb)
                    s2, o2 = decompose_index(i2, self.n_orb)
                    s3, o3 = decompose_index(i3, self.n_orb)
                    s4, o4 = decompose_index(i4, self.n_orb)
                else:
                    s1, o1 = 0, i1
                    s2, o2 = 0, i2
                    s3, o3 = 0, i3
                    s4, o4 = 0, i4

                s12 = block2.get_index(icrsh, s1, icrsh, s2)
                s43 = block2.get_index(icrsh, s4, icrsh, s3)
                inner12 = inner2.get_index(o1, o2)
                inner43 = inner2.get_index(o4, o3)

                if (s12, s43) not in x_loc:
                    x_loc[(s12, s43)] = numpy.zeros((n_inner**2, n_inner**2, 2*n_w2f, 2*n_w2f), dtype=complex)
                x_loc[(s12, s43)][inner12, inner43, :, :] = data_wb[:, :]

            # save
            h5bse.save(key=('X_loc', wb), data=x_loc)

    def name(self):
        return "pomerol"
