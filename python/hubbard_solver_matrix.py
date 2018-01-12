# Original file: https://github.com/TRIQS/hubbardI.git
#                python/hubbard_solver.py
# commit d919676f472bd5855751159bb91d1737ee9a46ea
# Date:   Tue Nov 21 19:19:26 2017 +0100

################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

from types import *
from pytriqs.operators.util.U_matrix import *
from pytriqs.gf.local import *
from pytriqs.applications.impurity_solvers.hubbard_I.hubbard_I import gf_hi_fullu, sigma_atomic_fullu
import pytriqs.utility.mpi as mpi
from itertools import izip
import numpy
import copy


class Solver:
    """
       Hubbard I Solver
    """

    # initialisation:
    def __init__(self, beta, norb, n_msb=1025, use_spin_orbit=False, nmoments=5):

        self.name = "Hubbard I"
        self.beta = beta
        self.Nmsb = n_msb
        self.UseSpinOrbit = use_spin_orbit
        self.Converged = False
        self.Nspin = 2
        self.Nmoments = nmoments
        self.Nlm = norb

        if use_spin_orbit:
            # no blocks!
            self.gf_struct = [('ud', range(2*self.Nlm))]
        else:
            # up/down blocks:
            self.gf_struct = [('up', range(self.Nlm)), ('down', range(self.Nlm))]

        # construct Greens functions:
        self.a_list = [a for a, al in self.gf_struct]

        def glist():
            return [GfImFreq(indices=ind, beta=self.beta, n_points=self.Nmsb) for blc, ind in self.gf_struct]
        self.G_iw = BlockGf(name_list=self.a_list, block_list=glist(), make_copies=False)
        self.G_Old = self.G_iw.copy()
        self.G0_iw = self.G_iw.copy()
        self.Sigma_iw = self.G_iw.copy()
        self.Sigma_Old = self.G_iw.copy()

        # prepare self.ealmat
        self.ealmat = numpy.zeros([self.Nlm*self.Nspin, self.Nlm*self.Nspin], numpy.complex_)

        # Define Atomic Levels Dictionary according to the GF Bloc Structure
        self.Eff_Atomic_Levels = {}
        for a, al in self.gf_struct:
            if self.UseSpinOrbit:
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm*2, self.Nlm*2], numpy.complex_)
            else:
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm, self.Nlm], numpy.complex_)

    def solve(self, u_mat, verbosity=0, iteration_number=1, test_convergence=0.0001,
              n_lev=0, remove_split=False):
        """Calculation of the impurity Greens function using Hubbard-I"""

        if self.Converged:
            mpi.report("Solver %(name)s has already converged: SKIPPING" % self.__dict__)
            return

        if mpi.is_master_node():
            my_verbosity = verbosity
        else:
            my_verbosity = 0

        # self.Nmoments = 5

        u_para, u_antipara = reduce_4index_to_2index(u_mat)

        mesh = [x for x in self.G_iw.mesh]
        zmsb = numpy.array([x for x in mesh], numpy.complex_)

        # for the tails:
        tailtempl = {}
        for sig, g in self.G_iw:
            tailtempl[sig] = copy.deepcopy(g.tail)
            for i in range(9):
                tailtempl[sig][i] *= 0.0

        self.__save_eal('eal.dat', iteration_number)

        mpi.report("Starting Fortran solver %(name)s" % self.__dict__)

        self.Sigma_Old <<= self.Sigma_iw
        self.G_Old <<= self.G_iw

        # call the fortran solver:
        temp = 1.0/self.beta
        gf, tail, atocc, atmag = gf_hi_fullu(e0f=self.ealmat, ur=u_mat, umn=u_para, ujmn=u_antipara,
                                             zmsb=zmsb, nmom=self.Nmoments, ns=self.Nspin, temp=temp,
                                             verbosity=my_verbosity, remove_split=remove_split,
                                             nlev_cf=n_lev)

        # self.sig = sigma_atomic_fullu(gf=self.gf, e0f=self.eal, zmsb=self.zmsb, ns=self.Nspin, nlm=self.Nlm)

        if my_verbosity == 0:
            # No fortran output, so give basic results here
            mpi.report("Atomic occupancy in Hubbard I Solver  : %s" % atocc)
            mpi.report("Atomic magn. mom. in Hubbard I Solver : %s" % atmag)

        # transfer the data to the GF class:
        if self.UseSpinOrbit:
            nlmtot = self.Nlm*2         # only one block in this case!
        else:
            nlmtot = self.Nlm

        g_trans = {}
        isp = -1
        for a, al in self.gf_struct:
            isp += 1
            g_trans[a] = numpy.array(
                gf[isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot, :]).transpose((2, 0, 1)).copy()
            for i in range(min(self.Nmoments, 8)):
                tailtempl[a][i+1] = tail[i][isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot]

        def glist():
            return [GfImFreq(indices=ind, beta=self.beta, n_points=self.Nmsb) for blc, ind in self.gf_struct]
        self.G_iw = BlockGf(name_list=self.a_list, block_list=glist(), make_copies=False)

        self.__copy_gf(self.G_iw, g_trans, tailtempl)

        # Self energy:
        self.G0_iw <<= iOmega_n

        eal0 = [self.ealmat[isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot] for isp in range((2*self.Nlm)/nlmtot)]
        self.G0_iw -= eal0
        self.Sigma_iw <<= self.G0_iw - inverse(self.G_iw)

        # invert G0
        self.G0_iw.invert()

        def test_distance(gf1, gf2, dist):
            def f(gf01, gf02):
                # print abs(G1.data - G2.data)
                dist_gf = max(abs(gf01.data - gf02.data).flatten())
                abs_gf = max(abs(gf01.data).flatten())
                return dist_gf <= abs_gf*dist
            return reduce(lambda x1, y1: x1 and y1, [f(g1, g2) for (i1, g1), (i2, g2) in izip(gf1, gf2)])

        mpi.report("\nChecking Sigma for convergence...\nUsing tolerance %s" % test_convergence)
        self.Converged = test_distance(self.Sigma_iw, self.Sigma_Old, test_convergence)

        if self.Converged:
            mpi.report("Solver HAS CONVERGED")
        else:
            mpi.report("Solver has not yet converged")

    def gf_realomega(self, ommin, ommax, n_om, u_mat, verbosity=0, broadening=0.01, n_lev=0, remove_split=False):
        """Calculates the GF and spectral function on the real axis."""

        delta_om = (ommax-ommin)/(1.0*(n_om-1))

        omega = numpy.zeros([n_om], numpy.complex_)

        u_para, u_antipara = reduce_4index_to_2index(u_mat)

        for i in range(n_om):
            omega[i] = ommin + delta_om * i + 1j * broadening

        tailtempl = {}
        for sig, g in self.G_iw:
            tailtempl[sig] = copy.deepcopy(g.tail)
            for i in range(9):
                tailtempl[sig][i] *= 0.0

        temp = 1.0/self.beta
        gf, tail, atocc, atmag = gf_hi_fullu(e0f=self.ealmat, ur=u_mat, umn=u_para, ujmn=u_antipara,
                                             zmsb=omega, nmom=self.Nmoments, ns=self.Nspin, temp=temp,
                                             verbosity=verbosity, remove_split=remove_split, nlev_cf=n_lev)

        # transfer the data to the GF class:
        if self.UseSpinOrbit:
            nlmtot = self.Nlm*2         # only one block in this case!
        else:
            nlmtot = self.Nlm

        g_trans = {}
        isp = -1
        for a, al in self.gf_struct:
            isp += 1
            g_trans[a] = numpy.array(
                gf[isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot, :]).transpose((2, 0, 1)).copy()
            for i in range(min(self.Nmoments, 8)):
                tailtempl[a][i+1] = tail[i][isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot]

        def glist():
            return [GfReFreq(indices=ind, window=(ommin, ommax), n_points=n_om) for blc, ind in self.gf_struct]
        self.G_iw = BlockGf(name_list=self.a_list, block_list=glist(), make_copies=False)

        self.__copy_gf(self.G_iw, g_trans, tailtempl)

        # Self energy:
        self.G0_iw = self.G_iw.copy()
        self.Sigma_iw = self.G_iw.copy()
        self.G0_iw <<= Omega + 1j*broadening

        eal0 = [self.ealmat[isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot] for isp in range((2*self.Nlm)/nlmtot)]
        self.G0_iw -= eal0
        self.Sigma_iw <<= self.G0_iw - inverse(self.G_iw)
        self.Sigma_iw.note = 'ReFreq'          # This is important for the put_Sigma routine!!!

        # sigmamat = sigma_atomic_fullu(gf = gf, e0f = self.ealmat, zmsb = omega, nlm = self.Nlm, ns = self.Nspin)

        # return omega, gf, sigmamat

    def __save_eal(self, filename, it):
        if mpi.is_master_node():
            f = open(filename, 'a')
            f.write('\neff. atomic levels, Iteration %s\n' % it)
            for i in range(self.Nlm*self.Nspin):
                for j in range(self.Nlm*self.Nspin):
                    f.write("%10.6f %10.6f   " % (self.ealmat[i, j].real, self.ealmat[i, j].imag))
                f.write("\n")
            f.close()

    def __copy_gf(self, gf, data, tail):
        """ Copies data and tail to Gf object GF """
        for s, g in gf:
            g.data[:, :, :] = data[s][:, :, :]
            for imom in range(1, min(self.Nmoments, 8)):
                g.tail.data[1+imom, :, :] = tail[s][imom]

    def set_atomic_levels(self, eal):
        """ Helps to set correctly the variables for the atomic levels from a dictionary."""

        assert (type(eal) == DictType), "Give a dictionary to set_atomic_levels!"

        cnt = 0
        self.ealmat[:, :] *= 0.0

        if self.UseSpinOrbit:
            self.Eff_Atomic_Levels["ud"] = copy.deepcopy(eal["ud"])

            for ii in range(self.Nlm*2):
                for jj in range(self.Nlm*2):
                    self.ealmat[ii, jj] = self.Eff_Atomic_Levels["ud"][ii, jj]

        else:
            for ind in ["up", "down"]:
                self.Eff_Atomic_Levels[ind] = copy.deepcopy(eal[ind])

                for ii in range(self.Nlm):
                    for jj in range(self.Nlm):
                        self.ealmat[cnt*self.Nlm + ii, cnt*self.Nlm + jj] = self.Eff_Atomic_Levels[ind][ii, jj]
                cnt += 1


#        for ind in eal:
#            self.Eff_Atomic_Levels[ind] = copy.deepcopy(eal[ind])
#
#            if self.UseSpinOrbit:
#                for ii in range(self.Nlm*2):
#                    for jj in range(self.Nlm*2):
#                        self.ealmat[ii, jj] = self.Eff_Atomic_Levels[ind][ii, jj]
#            else:
#                for ii in range(self.Nlm):
#                    for jj in range(self.Nlm):
#                        self.ealmat[cnt*self.Nlm + ii, cnt*self.Nlm + jj] = self.Eff_Atomic_Levels[ind][ii, jj]
#
#            cnt += 1
