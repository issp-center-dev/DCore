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
import sys
from collections import namedtuple
import shlex
import json
import ast

from dcore._dispatcher import *

from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0, extract_bath_params, expand_path
from .base import SolverBase
from .pomerol import assign_from_numpy_array


class ScipySolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.
        """

        super(ScipySolver, self).__init__(beta, gf_struct, u_mat, n_iw)

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
        # exec_path = expand_path(params_kw['exec_path'])
        exec_path = os.path.join(os.path.dirname(__file__), 'scipy_sparse_main.py')
        # check_version(mpirun_command, exec_path)

        # bath fitting
        n_bath = params_kw.get('n_bath', 0)  # 0 for Hubbard-I approximation
        ph_symmetric = params_kw.get('ph_symmetric', False)
        fit_params = {}
        for key in ['fit_gtol',]:
            if key in params_kw:
                fit_params[key] = params_kw[key]

        n_sites = self.n_orb + n_bath

        # parameters from input
        n_eigen = params_kw.get('n_eigen', 100)  # number of eigenstates to be computed
        dim_full_diag = params_kw.get('dim_full_diag', 10000)
        weight_threshold = params_kw.get('weight_threshold', 1e-6)
        ncv = params_kw.get('ncv', None)
        # np = params_kw.get('np', 1)  # MPI
        eigen_solver = params_kw.get('eigen_solver', 'eigsh')
        gf_solver = params_kw.get('gf_solver', 'bicgstab')
        # gf_rtol = params_kw.get('gf_rtol', 1e-5)
        # gf_atol = params_kw.get('gf_atol', 0.0)
        check_n_eigen = params_kw.get('check_n_eigen', True)
        check_orthonormality = params_kw.get('check_orthonormality', True)

        # convert particle_numbers to a list
        particle_numbers_in = params_kw.get('particle_numbers', 'auto')
        if particle_numbers_in in ('auto', 'all'):
            particle_numbers = particle_numbers_in
        else:
            try:
                particle_numbers = ast.literal_eval(particle_numbers_in)
                assert isinstance(particle_numbers, (list, tuple))
            except:
                print(f"Error in parsing particle_numbers={particle_numbers_in!r}", file=sys.stderr)
                print(f"Allowed inputs are 'auto', 'all', or '[1,2,3]' etc.", file=sys.stderr)
                sys.exit(1)

        # fixed parameters
        file_input = "input.in"
        file_h0 = "h0.in"
        file_umat = "umat.in"
        # file_gf = "gf.dat"

        params_solver = {
            'n_sites': n_sites,
            'n_flavors': int(self.n_flavors),  # avoid error in json.dump (*)
            'n_orb': self.n_orb,
            'beta': self.beta,
            'n_eigen': n_eigen,  # number of eigenstates to be computed
            'flag_spin_conserve': 1 if not self.use_spin_orbit else 0,
            'file_h0': file_h0 + '.npy',
            'file_umat': file_umat + '.npy',
            'flag_gf': 1,
            # 'file_gf': file_gf,
            'n_iw': self.n_iw,
            'n_bath': n_bath,
            'dim_full_diag': dim_full_diag,
            'particle_numbers': particle_numbers,
            'weight_threshold': weight_threshold,
            'ncv': ncv,
            'eigen_solver': eigen_solver,
            'gf_solver': gf_solver,
            # 'gf_rtol': gf_rtol,
            # 'gf_atol': gf_atol,
            'check_n_eigen' : check_n_eigen,
            'check_orthonormality' : check_orthonormality,
        }
        # (*) TypeError: Object of type int64 is not JSON serializable

        with open(file_input, "w") as f:
            json.dump(params_solver, f, indent=2)

        # -------------------------------------------------------------------------

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

        n_bath_so = n_bath
        if self.use_spin_orbit:
            n_bath_so *= 2  # for up and down spins

        bath_levels, bath_hyb = extract_bath_params(self._Delta_iw, self.beta, self.block_names, n_bath_so, ph_symmetric=ph_symmetric, **fit_params)
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

        # Save H0
        numpy.save(file_h0, h0_full)

        # (1c) Set U_{ijkl} for the solver

        # Save U_{ijkl}
        numpy.save(file_umat, self.u_mat)

        # -------------------------------------------------------------------------

        # (2) Run a working horse
        print("\nSolving the impurity problem...")
        with open('./stdout.log', 'w') as output_f:
            # launch_mpi_subprocesses(mpirun_command, [exec_path, file_input], output_f)
            launch_mpi_subprocesses("", ["python3", exec_path, file_input], output_f)

        print("\nFinish impurity problem")

        # -------------------------------------------------------------------------

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

        # Load Gf(iw)
        gf = np.load("gf.npy")
        assert gf.shape == (self.n_flavors, self.n_flavors, self.n_iw)
        assert gf.dtype == complex

        if not self.use_spin_orbit:
            gf = gf.reshape((2, self.n_orb, 2, self.n_orb, self.n_iw))
            gf = numpy.einsum("sisjw->sijw", gf)  # delete spin off-diagonal elements
            assert gf.shape == (2, self.n_orb, self.n_orb, self.n_iw)
        else:
            gf = gf.reshape((1, self.n_flavors, self.n_flavors, self.n_iw))
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

    def name(self):
        return "scipy/sparse"
