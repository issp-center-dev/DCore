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
from itertools import product
import os
import subprocess

from triqs.gf import *
from triqs.operators import *

from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0, extract_bath_params
from .base import SolverBase

VERSION_REQUIRED = 1.5


def check_version(_mpirun_command, _exec_path):
    print(" Checking version...")
    with open('./version', 'w') as output_f:
        launch_mpi_subprocesses(_mpirun_command, [_exec_path, "--version"], output_f)
    with open('./version', 'r') as output_f:
        out = output_f.read()
    # pomerol2dcore version 1.0
    print(" |", out, end="")
    print(" | version required", VERSION_REQUIRED)
    version = float(out.split()[2])
    if version >= VERSION_REQUIRED:
        print(" OK")
    else:
        print(" ERROR: requirement not satisfied")
        exit(1)


def assign_from_numpy_array(g_block, data, block_names):

    # DCore:
    #   g_block[].data.shape = (2*n_iw, n_orb, n_orb)
    #                          (2*n_iw, 2*n_orb, 2*n_orb)
    #       Both positive and negative freq are stored

    # pomerol:
    #   data.shape = (2, n_orb, n_orb, n_iw)       w/ spin-orbit
    #                (1, 2*n_orb, 2*n_orb, n_iw)   w/o spin-orbit
    #       Only positive freq

    for i, bname in enumerate(block_names):
        gf = g_block[bname]

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
    spn = index // n_orb
    orb = index % n_orb
    return spn, orb


class PomerolSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.
        """

        super(PomerolSolver, self).__init__(beta, gf_struct, u_mat, n_iw)

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
        #   self.use_spin_orbit

        # print("params_kw =", params_kw)
        exec_path = os.path.expandvars(params_kw['exec_path'])
        check_version(mpirun_command, exec_path)

        # bath fitting
        n_bath = params_kw.get('n_bath', 0)  # 0 for Hubbard-I approximation
        fit_params = {}
        for key in ['fit_gtol',]:
            if key in params_kw:
                fit_params[key] = params_kw[key]

        # for BSE
        flag_vx = params_kw.get('flag_vx', 0)
        n_w2f = params_kw.get('n_w2f', None)
        n_w2b = params_kw.get('n_w2b', None)
        file_freqs = params_kw.get('file_freqs', None)

        # dynamical susceptibility
        flag_suscep = params_kw.get('flag_suscep', 0)
        n_wb = params_kw.get('n_wb', None)

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
            'n_wf': self.n_iw,
            'flag_vx': flag_vx,
            'n_w2f': n_w2f,
            'n_w2b': n_w2b,
            'file_freqs': file_freqs,
            'flag_suscep': flag_suscep,
            'n_wb': n_wb,
            'n_bath': n_bath,
        }

        with open(file_pomerol, "w") as f:
            for key, val in list(params_pomerol.items()):
                if val is not None:
                    print(key, "=", val, file=f)

        # (1a) If H0 is necessary:
        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin1, ..., spin2, spin2, ...
        h0_mat = extract_H0(self._G0_iw, self.block_names)
        assert h0_mat.shape == (self.n_flavors, self.n_flavors)

        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)
        # H0 is extracted from the tail of the Green's function.
        self._Delta_iw = delta(self._G0_iw)

        bath_levels, bath_hyb = extract_bath_params(self._Delta_iw, self.beta, self.block_names, n_bath, **fit_params)
        assert bath_levels.shape == (2*n_bath,)
        assert bath_hyb.shape == (self.n_flavors, 2*n_bath)

        # Construct (impurity + bath) Hamiltonian matrix of size (L1+L2) times (L1+L2)
        L1 = self.n_flavors
        L2 = 2*n_bath
        h0_full = numpy.zeros((L1 + L2, L1 + L2), dtype=complex)
        h0_full[0:L1, 0:L1] = h0_mat
        h0_full[0:L1, L1:L1+L2] = bath_hyb
        h0_full[L1:L1+L2, 0:L1] = bath_hyb.conj().T
        h0_full[L1:L1+L2,L1:L1+L2] = numpy.diag(bath_levels)

        # save H0 into a file
        with open(file_h0, "w") as f:
            for i, j in product(list(range(h0_full.shape[0])), list(range(h0_full.shape[1]))):
                # TODO: real or complex
                if abs(h0_full[i, j]) != 0:
                    print(i, j, h0_full[i,j].real, h0_full[i,j].imag, file=f)
                    # print(i, j, h0_full[i, j].real, file=f)

        # (1c) Set U_{ijkl} for the solver
        with open(file_umat, "w") as f:
            for i, j, k, l in product(list(range(self.n_flavors)), repeat=4):
                # TODO: real or complex
                if abs(self.u_mat[i, j, k, l]) != 0:
                    print(i, j, k, l, self.u_mat[i, j, k, l].real, self.u_mat[i, j, k, l].imag, file=f)
                    #print(i, j, k, l, self.u_mat[i, j, k, l].real, file=f)

        # (2) Run a working horse
        with open('./stdout.log', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, file_pomerol], output_f)

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

        # load data as a complex type
        gf_1d = numpy.loadtxt(file_gf).view(numpy.complex128).reshape(-1)
        if not self.use_spin_orbit:
            gf = gf_1d.reshape((2, self.n_orb, self.n_orb, self.n_iw))
        else:
            gf = gf_1d.reshape((1, self.n_flavors, self.n_flavors, self.n_iw))
        assign_from_numpy_array(self._Gimp_iw, gf, self.block_names)

        # Compute Sigma_iw
        # NOTE:
        #   compute G0(iw) from h0_mat instead of using self._G0_iw, because
        #   self._G0_iw includes more information than that passed to the
        #   solver (see extract_H0 for details).

        # Rearrange indices (imp_up, imp_dn, bath_up, bath_dn) --> (imp_up, bath_up, imp_dn bath_dn)
        index_order =  list(range(self.n_orb))                                      # imp_up
        index_order += list(range(2*self.n_orb, 2*self.n_orb + n_bath))             # bath_up
        index_order += list(range(self.n_orb, 2*self.n_orb))                        # imp_dn
        index_order += list(range(2*self.n_orb + n_bath, 2*self.n_orb + 2*n_bath))  # bath_dn
        index_order = numpy.array(index_order)
        h0_updn = h0_full[index_order, :][:, index_order]

        # Cut H0 into block structure
        n_block = len(self.gf_struct)
        n_inner = h0_full.shape[0] // n_block
        h0_block = [h0_updn[s*n_inner:(s+1)*n_inner, s*n_inner:(s+1)*n_inner] for s in range(n_block)]

        # Construct G0 including bath sites
        bath_names = [ "bath" + str(i_bath) for i_bath in range(n_bath)]
        gf_struct_full = { block: list(inner_names) + bath_names for block, inner_names in list(self.gf_struct.items())}
        g0_full = make_block_gf(GfImFreq, gf_struct_full, self.beta, self.n_iw)
        g0_full << iOmega_n
        for i, block in enumerate(self.block_names):
            g0_full[block] -= h0_block[i]
        g0_full.invert()

        # Project G0 onto impurity site
        g0_imp = make_block_gf(GfImFreq, self.gf_struct, self.beta, self.n_iw)
        for block in self.block_names:
            for o1, o2 in product(self.gf_struct[block], repeat=2):
                g0_imp[block].data[:, o1, o2] = g0_full[block].data[:, o1, o2]

        self._Sigma_iw << inverse(g0_imp) - inverse(self._Gimp_iw)

    def _read_common(self, dir_name, fac=1.0):
        print("\n Reading files in directory '%s'" % dir_name)
        g2_loc = {}
        for i1, i2, i3, i4 in product(list(range(self.n_flavors)), repeat=4):
            filename = dir_name + "/%d_%d_%d_%d.dat" % (i1, i2, i3, i4)
            if not os.path.exists(filename):
                continue
            # print(" reading", filename)

            # load data as a complex type
            data = numpy.loadtxt(filename).view(numpy.complex128).reshape(-1)

            g2_loc[(i1, i2, i3, i4)] = data * fac
        print(" %d data loaded" % len(g2_loc))

        return g2_loc

    def _read_g2loc(self, params_kw):
        dir_g2 = params_kw.get('dir_g2', './two_particle')
        return self._read_common(dir_g2, 1./self.beta)

    def _read_chiloc(self, params_kw):
        dir_suscep = params_kw.get('dir_suscep', './susceptibility')
        return self._read_common(dir_suscep)

    def calc_Xloc_ph(self, rot, mpirun_command, num_wf, num_wb, params_kw):
        """
        compute local G2 in p-h channel
            X_loc = < c_{i1}^+ ; c_{i2} ; c_{i4}^+ ; c_{i3} >

        Parameters
        ----------
        rot
        mpirun_command
        num_wf
        num_wb
        params_kw

        Returns
        -------
        G2_loc : dict
            key = (i1, i2, i3, i4)
            val = numpy.ndarray(n_w2b, 2*n_w2f, 2*n_w2f)

        chi_loc : dict (None if not computed)
            key = (i1, i2, i3, i4)
            val = numpy.ndarray(n_w2b)
        """

        params_kw['flag_vx'] = 1
        params_kw['n_w2f'] = num_wf
        params_kw['n_w2b'] = num_wb
        params_kw['flag_suscep'] = 1
        params_kw['n_wb'] = num_wb

        self.solve(rot, mpirun_command, params_kw)

        g2_loc = self._read_g2loc(params_kw)
        # 1d array --> (wb, wf1, wf2)
        for key, data in list(g2_loc.items()):
            g2_loc[key] = data.reshape((num_wb, 2*num_wf, 2*num_wf))

        chi_loc = self._read_chiloc(params_kw)

        return g2_loc, chi_loc


    def calc_G2loc_ph_sparse(self, rot, mpirun_command, wsample_ph, params_kw):
        """

        Parameters
        ----------
        wsample_ph : 3*numpy.ndarray[N, 3]
           Sampling frequencies (in the order of fermion, fermion, boson)

        Returns
        -------
        g2_loc : complex array of shape (nflavor, nflavor, nflavor, nflavor, nfreq)
        """

        # Save frequencies list
        file_freqs = "freqs.in"
        with open(file_freqs, "w") as f:
            for wf1, wf2, wb in zip(*wsample_ph):
                print(wb//2, wf1//2, wf2//2, file=f)

        params_kw['flag_vx'] = 1
        params_kw['file_freqs'] = file_freqs
        params_kw['flag_suscep'] = 1

        self.solve(rot, mpirun_command, params_kw)

        G2loc_dict = self._read_g2loc(params_kw)
        G2loc = numpy.zeros((self.n_flavors,) * 4 + (wsample_ph[0].size,), dtype=numpy.complex128)
        # To irbasis_x's convention
        for k, v in G2loc_dict.items():
            G2loc[k[1], k[0], k[2], k[3], :] = self.beta * v
        return G2loc


    def name(self):
        return "pomerol"
