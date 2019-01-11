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
from pytriqs.gf.local import *
from pytriqs.applications.impurity_solvers.hubbard_I.hubbard_I import gf_hi_fullu, sigma_atomic_fullu
import pytriqs.utility.mpi as mpi
from itertools import izip
import numpy
import copy

class HubbardISolver(SolverBase):
    """
       Hubbard I Solver
    """

    # initialisation:
    def __init__(self, beta, gf_struct, u_mat, n_iw=1025, n_tau=10001):

        super(HubbardISolver, self).__init__(beta, gf_struct, u_mat, n_iw, n_tau)

        self.name = "Hubbard I"
        self.beta = beta
        self.n_iw = n_iw
        self.UseSpinOrbit = self.use_spin_orbit # From SolverBase
        self.Converged = False
        self.Nspin = 2
        self.Nmoments = 5
        self.Nlm = self.n_orb # From SolverBase

        if self.UseSpinOrbit:
            # no blocks!
            self.gf_struct = [('ud', range(2*self.Nlm))]
        else:
            # up/down blocks:
            self.gf_struct = [('up', range(self.Nlm)), ('down', range(self.Nlm))]

        # construct Greens functions:
        self.a_list = [a for a, al in self.gf_struct]

        def glist():
            return [GfImFreq(indices=ind, beta=self.beta, n_points=self.n_iw) for blc, ind in self.gf_struct]
        self.G_Old = self._Gimp_iw.copy()
        self.Sigma_Old = self._Gimp_iw.copy()

        # prepare self.ealmat
        self.ealmat = numpy.zeros([self.Nlm*self.Nspin, self.Nlm*self.Nspin], numpy.complex_)

        # Define Atomic Levels Dictionary according to the GF Bloc Structure
        self.Eff_Atomic_Levels = {}
        for a, al in self.gf_struct:
            if self.UseSpinOrbit:
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm*2, self.Nlm*2], numpy.complex_)
            else:
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm, self.Nlm], numpy.complex_)

    def solve(self, **params):
        #
        u_mat_spin_full = params['u_mat']


        u_mat = numpy.zeros((self.n_orb, self.n_orb, self.n_orb, self.n_orb), numpy.complex)
        u_mat[:, :, :, :] = u_mat_spin_full[0:self.n_orb, 0:self.n_orb, 0:self.n_orb, 0:self.n_orb]

        verbosity = params['verbosity'] if 'verbosity' in params else 0
        eal = {}
        for bname, gf in self._G0_iw:
            eal[bname] = numpy.array(gf.tail[2])
        self.set_atomic_levels(eal)

        test_convergence=0.0001
        n_lev=0
        remove_split=False

        if self.Converged:
            mpi.report("Solver %(name)s has already converged: SKIPPING" % self.__dict__)
            return

        if mpi.is_master_node():
            my_verbosity = verbosity
        else:
            my_verbosity = 0

        # self.Nmoments = 5

        u_para, u_antipara = reduce_4index_to_2index(u_mat)

        mesh = [x for x in self._Gimp_iw.mesh]
        zmsb = numpy.array([x for x in mesh], numpy.complex_)

        # for the tails:
        tailtempl = {}
        for sig, g in self._Gimp_iw:
            tailtempl[sig] = copy.deepcopy(g.tail)
            for i in range(9):
                tailtempl[sig][i] *= 0.0

        # self.__save_eal('eal.dat', iteration_number)

        print("Starting Fortran solver %(name)s" % self.__dict__)

        self.Sigma_Old <<= self._Sigma_iw
        self.G_Old <<= self._Gimp_iw

        # call the fortran solver:
        temp = 1.0/self.beta
        gf = numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm, len(zmsb)), numpy.complex_)
        tail = [numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm), numpy.complex_) for i in range(self.Nmoments)]
        atocc = 0.0
        atmag = 0.0
        gf, tail, atocc, atmag = gf_hi_fullu(e0f=self.ealmat, ur=u_mat, umn=u_para, ujmn=u_antipara,
                                             zmsb=zmsb, nmom=self.Nmoments, ns=self.Nspin, temp=temp,
                                             verbosity=my_verbosity, remove_split=remove_split,
                                             nlev_cf=n_lev)

        # self.sig = sigma_atomic_fullu(gf=self.gf, e0f=self.eal, zmsb=self.zmsb, ns=self.Nspin, nlm=self.Nlm)

        if my_verbosity == 0:
            # No fortran output, so give basic results here
            print("Atomic occupancy in Hubbard I Solver  : %s" % atocc)
            print("Atomic magn. mom. in Hubbard I Solver : %s" % atmag)

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
            return [GfImFreq(indices=ind, beta=self.beta, n_points=self.n_iw) for blc, ind in self.gf_struct]
        #self._Gimp_iw = BlockGf(name_list=self.a_list, block_list=glist(), make_copies=False)

        self.__copy_gf(self._Gimp_iw, g_trans, tailtempl)

        # Self energy:
        G0_iw = self._G0_iw.copy()
        G0_iw.zero()
        G0_iw <<= iOmega_n
        eal0 = [self.ealmat[isp*nlmtot:(isp+1)*nlmtot, isp*nlmtot:(isp+1)*nlmtot] for isp in range((2*self.Nlm)/nlmtot)]
        G0_iw -= eal0
        self._G0_iw << G0_iw
        self._Sigma_iw << self._G0_iw - inverse(self._Gimp_iw)

        # invert G0
        self._G0_iw.invert()

        def test_distance(gf1, gf2, dist):
            def f(gf01, gf02):
                # print abs(G1.data - G2.data)
                dist_gf = max(abs(gf01.data - gf02.data).flatten())
                abs_gf = max(abs(gf01.data).flatten())
                return dist_gf <= abs_gf*dist
            return reduce(lambda x1, y1: x1 and y1, [f(g1, g2) for (i1, g1), (i2, g2) in izip(gf1, gf2)])


    def gf_realomega(self, mesh, u_mat, verbosity=0, broadening=0.01, n_lev=0, remove_split=False):
        """Calculates the GF and spectral function on the real axis."""

        ommin, ommax, n_om = mesh

        delta_om = (ommax-ommin)/(1.0*(n_om-1))

        omega = numpy.zeros([n_om], numpy.complex_)

        u_para, u_antipara = reduce_4index_to_2index(u_mat)

        for i in range(n_om):
            omega[i] = ommin + delta_om * i + 1j * broadening

        tailtempl = {}
        for sig, g in self._Gimp_iw:
            tailtempl[sig] = copy.deepcopy(g.tail)
            for i in range(9):
                tailtempl[sig][i] *= 0.0

        temp = 1.0/self.beta
        gf = numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm, len(omega)), numpy.complex_)
        tail = [numpy.zeros((self.Nspin*self.Nlm, self.Nspin*self.Nlm), numpy.complex_) for i in range(self.Nmoments)]
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


    def __save_eal(self, filename, it):
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

    h_int = make_h_int(u_mat, gf_struct)

    # Create a working horse
    S = HubbardISolver(beta, gf_struct, u_mat, n_iw, n_tau)
    S.set_G0_iw(G0_iw)

    for k in params:
        # e.g. numpy.bool_ to bool
        params[k] = convert_to_built_in_scalar_type(params[k])

    params['u_mat'] =  u_mat
    S.solve(h_int=h_int, **params)

    if params['calc_Sigma_w']:
        S.gf_realomega(params['omeags'], u_mat)

    # Save results
    if mpi.is_master_node():
        with HDFArchive(os.path.abspath(output_file), 'w') as h:
            h['Sigma_iw'] = S.get_Sigma_iw()
            h['Gimp_iw'] = S.get_Gimp_iw()
            if params['calc_Sigma_w']:
                h['Sigma_w'] = S.get_Sigma_w()



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
