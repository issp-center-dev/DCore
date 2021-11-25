
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
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
##########################################################################

import copy
import numpy
from types import *
from ..gf import *
from ..h5 import *
from .. import mpi


class Symmetry:
    """
    This class provides the routines for applying symmetry operations for the k sums.
    It contains the permutations of the atoms in the unit cell, and the corresponding
    rotational matrices for each symmetry operation.
    """

    def __init__(self, hdf_file, subgroup=None):
        """
        Initialises the class.

        Parameters
        ----------
        hdf_file : string
                   Base name of the hdf5 archive with the symmetry data.
        subgroup : string, optional
                   Name of subgroup storing correlated-shell symmetry data. If not given, it is assumed that
                   the data is stored at the root of the hdf5 archive.
        """

        assert isinstance(hdf_file, str), "Symmetry: hdf_file must be a filename."
        self.hdf_file = hdf_file
        things_to_read = ['n_symm', 'n_atoms', 'perm',
                          'orbits', 'SO', 'SP', 'time_inv', 'mat', 'mat_tinv']
        for it in things_to_read:
            setattr(self, it, 0)

        if mpi.is_master_node():
            # Read the stuff on master:
            with HDFArchive(hdf_file, 'r') as ar:
                if subgroup is None:
                    ar2 = ar
                else:
                    ar2 = ar[subgroup]

                for it in things_to_read:
                    setattr(self, it, ar2[it])
            del ar2

        # Broadcasting
        for it in things_to_read:
            setattr(self, it, mpi.bcast(getattr(self, it)))

        # now define the mapping of orbitals:
        # self.orb_map[iorb] = jorb gives the permutation of the orbitals as given in the list, when the
        # permutation of the atoms is done:
        self.n_orbits = len(self.orbits)
        self.orb_map = [[0 for iorb in range(
            self.n_orbits)] for i_symm in range(self.n_symm)]
        for i_symm in range(self.n_symm):
            for iorb in range(self.n_orbits):
                srch = copy.deepcopy(self.orbits[iorb])
                srch['atom'] = self.perm[i_symm][self.orbits[iorb]['atom'] - 1]
                self.orb_map[i_symm][iorb] = self.orbits.index(srch)

    def symmetrize(self, obj):
        """
        Symmetrizes a given object. 

        Parameters
        ----------
        obj : list
              object to symmetrize. It has to be given as list, where its length is determined by the number 
              of equivalent members of the object. Two types of objects are supported:

              - BlockGf : list of Green's functions,
              - Matrices : The format is taken from density matrices as obtained from Green's functions (DictType).

        Returns
        -------
        symm_obj : list
                   Symmetrized object, of the same type as input object.
        """

        assert isinstance(
            obj, list), "symmetrize: obj has to be a list of objects."
        assert len(
            obj) == self.n_orbits, "symmetrize: obj has to be a list of the same length as defined in the init."

        if isinstance(obj[0], BlockGf):
            # here the result is stored, it is a BlockGf!
            symm_obj = [obj[i].copy() for i in range(len(obj))]
            for iorb in range(self.n_orbits):
                symm_obj[iorb].zero()      # set to zero
        else:
            # if not a BlockGf, we assume it is a matrix (density matrix), has
            # to be complex since self.mat is complex!
            symm_obj = [copy.deepcopy(obj[i]) for i in range(len(obj))]
            for iorb in range(self.n_orbits):
                if isinstance(symm_obj[iorb], dict):
                    for ii in symm_obj[iorb]:
                        symm_obj[iorb][ii] *= 0.0
                else:
                    symm_obj[iorb] *= 0.0

        for i_symm in range(self.n_symm):
            for iorb in range(self.n_orbits):
                l = self.orbits[iorb]['l']         # s, p, d, or f
                dim = self.orbits[iorb]['dim']
                jorb = self.orb_map[i_symm][iorb]

                if isinstance(obj[0], BlockGf):

                    tmp = obj[iorb].copy()
                    if self.time_inv[i_symm]:
                        tmp << tmp.transpose()
                    for bname, gf in tmp:
                        tmp[bname].from_L_G_R(self.mat[i_symm][iorb], tmp[bname], self.mat[
                                              i_symm][iorb].conjugate().transpose())
                    tmp *= 1.0 / self.n_symm
                    symm_obj[jorb] += tmp

                else:

                    if isinstance(obj[iorb], dict):
                        for ii in obj[iorb]:
                            if self.time_inv[i_symm] == 0:
                                symm_obj[jorb][ii] += numpy.dot(numpy.dot(self.mat[i_symm][iorb], obj[iorb][ii]),
                                                                self.mat[i_symm][iorb].conjugate().transpose()) / self.n_symm
                            else:
                                symm_obj[jorb][ii] += numpy.dot(numpy.dot(self.mat[i_symm][iorb], obj[iorb][ii].conjugate()),
                                                                self.mat[i_symm][iorb].conjugate().transpose()) / self.n_symm
                    else:
                        if self.time_inv[i_symm] == 0:
                            symm_obj[jorb] += numpy.dot(numpy.dot(self.mat[i_symm][iorb], obj[iorb]),
                                                        self.mat[i_symm][iorb].conjugate().transpose()) / self.n_symm
                        else:
                            symm_obj[jorb] += numpy.dot(numpy.dot(self.mat[i_symm][iorb], obj[iorb].conjugate()),
                                                        self.mat[i_symm][iorb].conjugate().transpose()) / self.n_symm

# Markus: This does not what it is supposed to do, check how this should work (keep for now)
#        if (self.SO == 0) and (self.SP == 0):
#            # add time inv:
            #mpi.report("Add time inversion")
#            for iorb in range(self.n_orbits):
#                if (isinstance(symm_obj[0],BlockGf)):
#                    tmp = symm_obj[iorb].copy()
#                    tmp << tmp.transpose()
#                    for bname,gf in tmp: tmp[bname].from_L_G_R(self.mat_tinv[iorb],tmp[bname],self.mat_tinv[iorb].transpose().conjugate())
#                    symm_obj[iorb] += tmp
#                    symm_obj[iorb] /= 2.0
#
#                else:
#                    if isinstance(symm_obj[iorb], dict):
#                        for ii in symm_obj[iorb]:
#                            symm_obj[iorb][ii] += numpy.dot(numpy.dot(self.mat_tinv[iorb],symm_obj[iorb][ii].conjugate()),
#                                                            self.mat_tinv[iorb].transpose().conjugate())
#                            symm_obj[iorb][ii] /= 2.0
#                    else:
#                        symm_obj[iorb] += numpy.dot(numpy.dot(self.mat_tinv[iorb],symm_obj[iorb].conjugate()),
#                                                    self.mat_tinv[iorb].transpose().conjugate())
#                        symm_obj[iorb] /= 2.0

        return symm_obj
