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

import argparse
import copy
import sys
import os
import numpy

from pytriqs.archive.hdf_archive import HDFArchive
import pytriqs.utility.mpi as mpi

from .base import make_h_int, SolverBase
from ..tools import convert_to_built_in_scalar_type

from types import *
from pytriqs.operators.util.U_matrix import *
from pytriqs.gf import *
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

    def solve(self, u_mat, verbosity=0, test_convergence=0.0001, n_lev=0, remove_split=False):
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

        # self.__save_eal('eal.dat', iteration_number)

        mpi.report("Starting Fortran solver %(name)s" % self.__dict__)

        self.Sigma_Old <<= self.Sigma_iw
        self.G_Old <<= self.G_iw

        # call the fortran solver:
        temp = 1.0/self.beta
        gf = numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm, len(zmsb)), numpy.complex_)
        tail = [numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm), numpy.complex_) for i in range(self.Nmoments)]
        atocc = 0.0
        atmag = 0.0
        if mpi.is_master_node():
            gf, tail, atocc, atmag = gf_hi_fullu(e0f=self.ealmat, ur=u_mat, umn=u_para, ujmn=u_antipara,
                                                 zmsb=zmsb, nmom=self.Nmoments, ns=self.Nspin, temp=temp,
                                                 verbosity=my_verbosity, remove_split=remove_split,
                                                 nlev_cf=n_lev)
        gf = mpi.bcast(gf)
        tail = mpi.bcast(tail)
        atocc = mpi.bcast(atocc)
        atmag = mpi.bcast(atmag)

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
        gf = numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm, len(omega)), numpy.complex_)
        tail = [numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm), numpy.complex_) for i in range(self.Nmoments)]
        if mpi.is_master_node():
            gf, tail, atocc, atmag = gf_hi_fullu(e0f=self.ealmat, ur=u_mat, umn=u_para, ujmn=u_antipara,
                                                 zmsb=omega, nmom=self.Nmoments, ns=self.Nspin, temp=temp,
                                                 verbosity=verbosity, remove_split=remove_split, nlev_cf=n_lev)
        gf = mpi.bcast(gf)
        tail = mpi.bcast(tail)

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
        self.G_w = BlockGf(name_list=self.a_list, block_list=glist(), make_copies=False)

        self.__copy_gf(self.G_w, g_trans, tailtempl)

        # Self energy:
        self.G0_w = self.G_w.copy()
        self.Sigma_w = self.G_w.copy()
        self.G0_w <<= Omega + 1j*broadening

        eal0 = [self.ealmat[isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot] for isp in range((2*self.Nlm)/nlmtot)]
        self.G0_w -= eal0
        self.Sigma_w <<= self.G0_w - inverse(self.G_w)
        self.Sigma_w.note = 'ReFreq'          # This is important for the put_Sigma routine!!!

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

def main(input_file, output_file):
    """
    Solve the impurity problem.
    """

    with HDFArchive(os.path.abspath(input_file), 'r') as h:
        rot = h['rot'] if 'rot' in h else None
        beta = h['beta']
        gf_struct = h['gf_struct']
        u_mat = h['u_mat']
        n_iw = h['n_iw']
        n_tau = h['n_tau']
        G0_iw = h['G0_iw']
        params = h['params']

    use_spin_orbit = len(gf_struct) == 1

    # Create a working horse
    norb = int(u_mat.shape[0]/2)
    S = Solver(beta=beta, norb=norb, n_msb=n_iw, use_spin_orbit=use_spin_orbit)
    eal = {}
    for bname, gf in G0_iw:
        eal[bname] = numpy.array(gf.tail[2])
    S.set_atomic_levels(eal)
    umat2 = u_mat[0:norb, 0:norb, 0:norb, 0:norb]
    #print("debug", umat2)
    #print("eal", eal)
    S.solve(u_mat=numpy.real(umat2), verbosity=True)

    calc_Sigma_w = 'calc_Sigma_w' in params and params['calc_Sigma_w']

    if calc_Sigma_w:
        mesh = (params['omega_min'], params['omega_max'], params['n_omega'])
        S.gf_realomega(ommin=params['omega_min'], ommax=params['omega_max'], n_om=params['n_omega'],
                              u_mat=numpy.real(umat2))

    # Save results
    if mpi.is_master_node():
        with HDFArchive(os.path.abspath(output_file), 'w') as h:
            h['Sigma_iw'] = S.Sigma_iw
            h['Gimp_iw'] = S.G_iw
            if calc_Sigma_w:
                h['Sigma_w'] = S.Sigma_w

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(description='Wrapper program for TRIQS/hubbard-I')
        parser.add_argument('input_file')
        parser.add_argument('output_file')
        args = parser.parse_args()
        main(args.input_file, args.output_file)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print("Unexpected error:", e)
        sys.exit(1)
