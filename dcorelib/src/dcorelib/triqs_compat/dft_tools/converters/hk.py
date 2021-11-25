
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn
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

from types import *
import numpy
from ...h5 import *
from ... import mpi
from math import sqrt
from .converter_tools import *


class HkConverter(ConverterTools):
    """
    Conversion from general H(k) file to an hdf5 file that can be used as input for the SumKDFT class.
    """

    def __init__(self, filename, hdf_filename=None, dft_subgrp='dft_input', symmcorr_subgrp='dft_symmcorr_input', repacking=False):
        """
        Initialise the class.

        Parameters
        ----------
        filename : string
                   Name of file containing the H(k) and other relevant data.
        hdf_filename : string, optional
                       Name of hdf5 archive to be created.
        dft_subgrp : string, optional
                     Name of subgroup storing necessary DFT data.
        symmcorr_subgrp : string, optional
                          Name of subgroup storing correlated-shell symmetry data.
                          The group is actually empty; it is just included for compatibility.
        repacking : boolean, optional
                    Does the hdf5 archive need to be repacked to save space?

        """

        assert isinstance(filename, str), "HkConverter: filename must be a filename."
        if hdf_filename is None:
            hdf_filename = filename + '.h5'
        self.hdf_file = hdf_filename
        self.dft_file = filename
        self.dft_subgrp = dft_subgrp
        self.symmcorr_subgrp = symmcorr_subgrp
        self.fortran_to_replace = {'D': 'E', '(': ' ', ')': ' ', ',': ' '}

        # Checks if h5 file is there and repacks it if wanted:
        import os.path
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)

    def convert_dft_input(self, first_real_part_matrix=True, only_upper_triangle=False, weights_in_file=False):
        """
        Reads the appropriate files and stores the data for the dft_subgrp in the hdf5 archive. 

        Parameters
        ----------
        first_real_part_matrix : boolean, optional
                                 Should all the real components for given k be read in first, followed by the imaginary parts?
        only_upper_triangle : boolean, optional
                              Should only the upper triangular part of H(k) be read in? 
        weights_in_file : boolean, optional
                          Are the k-point weights to be read in?

        """

        # Read and write only on the master node
        if not (mpi.is_master_node()):
            return
        mpi.report("Reading input from %s..." % self.dft_file)

        # R is a generator : each R.Next() will return the next number in the
        # file
        R = ConverterTools.read_fortran_file(
            self, self.dft_file, self.fortran_to_replace)
        try:
            # the energy conversion factor is 1.0, we assume eV in files
            energy_unit = 1.0
            # read the number of k points
            n_k = int(next(R))
            k_dep_projection = 0
            SP = 0                                        # no spin-polarision
            SO = 0                                        # no spin-orbit
            # total charge below energy window is set to 0
            charge_below = 0.0
            # density required, for setting the chemical potential
            density_required = next(R)
            symm_op = 0                                   # No symmetry groups for the k-sum

            # the information on the non-correlated shells is needed for
            # defining dimension of matrices:
            # number of shells considered in the Wanniers
            n_shells = int(next(R))
            # corresponds to index R in formulas
            # now read the information about the shells (atom, sort, l, dim):
            shell_entries = ['atom', 'sort', 'l', 'dim']
            shells = [{name: int(val) for name, val in zip(
                shell_entries, R)} for ish in range(n_shells)]

            # number of corr. shells (e.g. Fe d, Ce f) in the unit cell,
            n_corr_shells = int(next(R))
            # corresponds to index R in formulas
            # now read the information about the shells (atom, sort, l, dim, SO
            # flag, irep):
            corr_shell_entries = ['atom', 'sort', 'l', 'dim','SO','irep']
            corr_shells = [{name: int(val) for name, val in zip(
                corr_shell_entries, R)} for icrsh in range(n_corr_shells)]

            # determine the number of inequivalent correlated shells and maps,
            # needed for further reading
            [n_inequiv_shells, corr_to_inequiv,
                inequiv_to_corr] = ConverterTools.det_shell_equivalence(self, corr_shells)

            use_rotations = 0
            rot_mat = [numpy.identity(
                corr_shells[icrsh]['dim'], numpy.complex_) for icrsh in range(n_corr_shells)]
            rot_mat_time_inv = [0 for i in range(n_corr_shells)]

            # Representative representations are read from file
            n_reps = [1 for i in range(n_inequiv_shells)]
            dim_reps = [0 for i in range(n_inequiv_shells)]
            T = []
            for ish in range(n_inequiv_shells):
                # number of representatives ("subsets"), e.g. t2g and eg
                n_reps[ish] = int(next(R))
                dim_reps[ish] = [int(next(R)) for i in range(
                    n_reps[ish])]   # dimensions of the subsets

                # The transformation matrix:
                # is of dimension 2l+1, it is taken to be standard d (as in
                # Wien2k)
                ll = 2 * corr_shells[inequiv_to_corr[ish]]['l'] + 1
                lmax = ll * (corr_shells[inequiv_to_corr[ish]]['SO'] + 1)
                T.append(numpy.zeros([lmax, lmax], numpy.complex_))

                T[ish] = numpy.array([[0.0, 0.0, 1.0, 0.0, 0.0],
                                      [1.0 / sqrt(2.0), 0.0, 0.0,
                                       0.0, 1.0 / sqrt(2.0)],
                                      [-1.0 / sqrt(2.0), 0.0, 0.0,
                                       0.0, 1.0 / sqrt(2.0)],
                                      [0.0, 1.0 /
                                          sqrt(2.0), 0.0, -1.0 / sqrt(2.0), 0.0],
                                      [0.0, 1.0 / sqrt(2.0), 0.0, 1.0 / sqrt(2.0), 0.0]])

            # Spin blocks to be read:
            # number of spins to read for Norbs and Ham, NOT Projectors
            n_spin_blocs = SP + 1 - SO

            # define the number of n_orbitals for all k points: it is the
            # number of total bands and independent of k!
            n_orbitals = numpy.ones(
                [n_k, n_spin_blocs], numpy.int) * sum([sh['dim'] for sh in shells])

            # Initialise the projectors:
            proj_mat = numpy.zeros([n_k, n_spin_blocs, n_corr_shells, max(
                [crsh['dim'] for crsh in corr_shells]), numpy.max(n_orbitals)], numpy.complex_)

            # Read the projectors from the file:
            for ik in range(n_k):
                for icrsh in range(n_corr_shells):
                    for isp in range(n_spin_blocs):

                        # calculate the offset:
                        offset = 0
                        n_orb = 0
                        for ish in range(n_shells):
                            if (n_orb == 0):
                                if (shells[ish]['atom'] == corr_shells[icrsh]['atom']) and (shells[ish]['sort'] == corr_shells[icrsh]['sort']):
                                    n_orb = corr_shells[icrsh]['dim']
                                else:
                                    offset += shells[ish]['dim']

                        proj_mat[ik, isp, icrsh, 0:n_orb,
                                 offset:offset + n_orb] = numpy.identity(n_orb)

            # now define the arrays for weights and hopping ...
            # w(k_index),  default normalisation
            bz_weights = numpy.ones([n_k], numpy.float_) / float(n_k)
            hopping = numpy.zeros([n_k, n_spin_blocs, numpy.max(
                n_orbitals), numpy.max(n_orbitals)], numpy.complex_)

            if (weights_in_file):
                # weights in the file
                for ik in range(n_k):
                    bz_weights[ik] = next(R)

            # if the sum over spins is in the weights, take it out again!!
            sm = sum(bz_weights)
            bz_weights[:] /= sm

            # Grab the H
            for isp in range(n_spin_blocs):
                for ik in range(n_k):
                    n_orb = n_orbitals[ik, isp]

                    # first read all real components for given k, then read
                    # imaginary parts
                    if (first_real_part_matrix):

                        for i in range(n_orb):
                            if (only_upper_triangle):
                                istart = i
                            else:
                                istart = 0
                            for j in range(istart, n_orb):
                                hopping[ik, isp, i, j] = next(R)

                        for i in range(n_orb):
                            if (only_upper_triangle):
                                istart = i
                            else:
                                istart = 0
                            for j in range(istart, n_orb):
                                hopping[ik, isp, i, j] += next(R) * 1j
                                if ((only_upper_triangle)and(i != j)):
                                    hopping[ik, isp, j, i] = hopping[
                                        ik, isp, i, j].conjugate()

                    else:  # read (real,im) tuple

                        for i in range(n_orb):
                            if (only_upper_triangle):
                                istart = i
                            else:
                                istart = 0
                            for j in range(istart, n_orb):
                                hopping[ik, isp, i, j] = next(R)
                                hopping[ik, isp, i, j] += next(R) * 1j

                                if ((only_upper_triangle)and(i != j)):
                                    hopping[ik, isp, j, i] = hopping[
                                        ik, isp, i, j].conjugate()
            # keep some things that we need for reading parproj:
            things_to_set = ['n_shells', 'shells', 'n_corr_shells', 'corr_shells',
                             'n_spin_blocs', 'n_orbitals', 'n_k', 'SO', 'SP', 'energy_unit']
            for it in things_to_set:
                setattr(self, it, locals()[it])
        except StopIteration:  # a more explicit error if the file is corrupted.
            raise "HK Converter : reading file dft_file failed!"

        R.close()

        # Save to the HDF5:
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.dft_subgrp in ar):
                ar.create_group(self.dft_subgrp)
            things_to_save = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                          'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                          'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights', 'hopping',
                          'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
            for it in things_to_save:
                ar[self.dft_subgrp][it] = locals()[it]
