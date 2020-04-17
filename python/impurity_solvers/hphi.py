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

from ..pytriqs_gf_compat import *
# from pytriqs.archive import HDFArchive
from pytriqs.operators import *

from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0, extract_bath_params
from .base import SolverBase



class HPhiSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.
        """

        super(HPhiSolver, self).__init__(beta, gf_struct, u_mat, n_iw)

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

        exec_path = os.path.expandvars(params_kw['exec_path'])

        # bath fitting
        n_bath = params_kw.get('n_bath', 0)  # 0 for Hubbard-I approximation
        fit_params = {}
        for key in ['fit_gtol',]:
            if key in params_kw:
                fit_params[key] = params_kw[key]

        n_site = self.n_orb + n_bath

# -------------------------------------------------------------------------

        with open('./namelist.def', 'w') as f:
            print("""\
ModPara  modpara.def
CalcMod  calcmod.def
LocSpin  locspn.def
Trans  trans.def
InterAll  interall.def
            """, file=f)

# -------------------------------------------------------------------------

        with open('./modpara.def', 'w') as f:
            print("""\
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          {0}
Lanczos_max    2000
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     4.000000000000000e+00
NumAve         5
ExpecInterval  20
NOmega         200
OmegaMax       1.200000000000000e+01     4.000000000000000e-02
OmegaMin       -1.200000000000000e+01    4.000000000000000e-02
OmegaOrg       0.000000000000000e+00     0.000000000000000e+00\
""".format(n_site), file=f)

# -------------------------------------------------------------------------

        with open('./calcmod.def', 'w') as f:
            print("""\
#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution
#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC
#Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart
#CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save
CalcType   3
CalcModel   3
ReStart   0
CalcSpec   0
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   0
InputHam   0
OutputHam   0
OutputExVec   0\
""", file=f)

# -------------------------------------------------------------------------

        with open('./locspn.def', 'w') as f:
            print("""\
================================
NlocalSpin     0
================================
========i_0LocSpn_1IteElc ======
================================\
""", file=f)
            for i in range(n_site):
                print("{0} 0".format(i), file=f)

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

        bath_levels, bath_hyb = extract_bath_params(self._Delta_iw, self.beta, self.block_names, n_bath, **fit_params)
        assert bath_levels.shape == (2*n_bath,)
        assert bath_hyb.shape == (self.n_flavors, 2*n_bath)

        # number of hopping terms
        n_trans  = (2*self.n_orb)**2  # correlated sites
        n_trans += (2*n_bath)  # bath levels
        n_trans += 2*(2*self.n_orb)*(2*n_bath)  # hopping between correlated sites and bath sites

        with open('./trans.def', 'w') as f:
            print("""\
========================
NTransfer      {0}
========================
========i_j_s_tijs======
========================
""".format(n_trans), file=f)

            # correlated sites
            h0_isjs = h0_mat.reshape((2, self.n_orb, 2, self.n_orb))
            for s1, s2 in product(range(2), repeat=2):
                for i1, i2 in product(range(self.n_orb), repeat=2):
                    t = -h0_isjs[s1, i1, s2, i2]
                    print(i1, s1, i2, s2, t.real, t.imag, file=f)

            # bath levels
            bath_levels_is = bath_levels.reshape((2, n_bath))
            for s1 in range(2):
                for i1 in range(n_bath):
                    j1 = self.n_orb + i1
                    eps = -bath_levels_is[s1, i1]
                    print(j1, s1, j1, s1, eps.real, eps.imag, file=f)

            # hopping between correlated sites and bath sites
            bath_hyb_is = bath_hyb.reshape((2, self.n_orb, 2, n_bath))
            for s1 in range(2):
                for i1 in range(self.n_orb):
                    for s2 in range(2):
                        for i2 in range(n_bath):
                            j1 = i1
                            j2 = self.n_orb + i2
                            v = -bath_hyb_is[s1, i1, s2, i2]
                            print(j1, s1, j2, s2, v.real, v.imag, file=f)
                            print(j2, s2, j1, s1, v.real, -v.imag, file=f)

        # -------------------------------------------------------------------------

        # (1c) Set U_{ijkl} for the solver
        # for i, j, k, l in product(range(self.n_flavors), repeat=4):
        #     self.u_mat[i, j, k, l]

        n_inter = (2*self.n_orb)**4

        with open('./interall.def', 'w') as f:
            print("""\
======================
NInterAll      {0}
======================
========zInterAll=====
======================
""".format(n_inter), file=f)

        # (1/2) U_{1234} c_1^+ c_2^+ c_4 c_3
        # = I_{1324} c_1^+ c_3 c_2^+ c_4

            u_1234 = self.u_mat.reshape((2, self.n_orb, 2, self.n_orb, 2, self.n_orb, 2, self.n_orb))
            for s1, s2, s3, s4 in product(range(2), repeat=4):
                for o1, o2, o3, o4 in product(range(self.n_orb), repeat=4):
                    u = u_1234[s1, o1, s2, o2, s3, o3, s4, o4] / 2.
                    if s1==s2==s3==s4 and o1==o2==o3==o4:
                        u = 0.0
                    print(o1, s1, o3, s3, o2, s2, o4, s4, u.real, u.imag, file=f)

        # (2) Run a working horse
        with open('./stdout.log', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, '-e', 'namelist.def'], output_f)

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

    def name(self):
        return "HPhi"
