
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

###
#  Wannier90 to HDF5 converter for the SumkDFT class of dfttools/TRIQS;
#
#   written by Gabriele Sclauzero (Materials Theory, ETH Zurich), Dec 2015 -- Jan 2016,
#   under the supervision of Claude Ederer (Materials Theory).
#   Partially based on previous work by K. Dymkovski and the DFT_tools/TRIQS team.
#
#  Limitations of the current implementation:
# - the case with SO=1 is not considered at the moment
# - the T rotation matrices are not used in this implementation
# - projectors for uncorrelated shells (proj_mat_all) cannot be set
#
#  Things to be improved/checked:
# - the case with SP=1 might work, but was never tested (do we need to define
#   rot_mat_time_inv also if symm_op = 0?)
# - the calculation of rot_mat in find_rot_mat() relies on the eigenvalues of H(0);
#   this might fail in presence of degenerate eigenvalues (now just prints warning)
# - the FFT is always done in serial mode (because all converters run serially);
#   this can become very slow with a large number of R-vectors/k-points
# - make the code more MPI safe (error handling): if we run with more than one process
#   and an error occurs on the masternode, the calculation does not abort
###


from types import *
import numpy
import math
from ...h5 import *
from .converter_tools import *
from itertools import product
import os.path


class Wannier90Converter(ConverterTools):
    """
    Conversion from Wannier90 output to an hdf5 file that can be used as input for the SumkDFT class.
    """

    def __init__(self, seedname, hdf_filename=None, dft_subgrp='dft_input',
                 symmcorr_subgrp='dft_symmcorr_input', repacking=False):
        """
        Initialise the class.

        Parameters
        ----------
        seedname : string
            Base name of Wannier90 files
        hdf_filename : string, optional
            Name of hdf5 archive to be created
        dft_subgrp : string, optional
            Name of subgroup storing necessary DFT data
        symmcorr_subgrp : string, optional
            Name of subgroup storing correlated-shell symmetry data
        repacking : boolean, optional
            Does the hdf5 archive need to be repacked to save space?

        """

        self._name = "Wannier90Converter"
        assert isinstance(seedname, str), self._name + \
            ": Please provide the DFT files' base name as a string."
        if hdf_filename is None:
            hdf_filename = seedname + '.h5'
        self.hdf_file = hdf_filename
        # if the w90 output is seedname_hr.dat, the input file for the
        # converter must be called seedname.inp
        self.inp_file = seedname + '.inp'
        self.w90_seed = seedname
        self.dft_subgrp = dft_subgrp
        self.symmcorr_subgrp = symmcorr_subgrp
        self.fortran_to_replace = {'D': 'E'}
        # threshold below which matrix elements from wannier90 should be
        # considered equal
        self._w90zero = 2.e-6

        # Checks if h5 file is there and repacks it if wanted:
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)

    def convert_dft_input(self):
        """
        Reads the appropriate files and stores the data for the

        - dft_subgrp
        - symmcorr_subgrp

        in the hdf5 archive. 

        """

        # Read and write only on the master node
        if not (mpi.is_master_node()):
            return
        mpi.report("Reading input from %s..." % self.inp_file)

        # R is a generator : each R.Next() will return the next number in the
        # file
        R = ConverterTools.read_fortran_file(
            self, self.inp_file, self.fortran_to_replace)
        shell_entries = ['atom', 'sort', 'l', 'dim']
        corr_shell_entries = ['atom', 'sort', 'l', 'dim', 'SO', 'irep']
        # First, let's read the input file with the parameters needed for the
        # conversion
        try:
            # read k - point mesh generation option
            kmesh_mode = int(next(R))
            if kmesh_mode >= 0:
                # read k-point mesh size from input
                nki = [int(next(R)) for idir in range(3)]
            else:
                # some default grid, if everything else fails...
                nki = [8, 8, 8]
            # read the total number of electrons per cell
            density_required = float(next(R))
            # we do not read shells, because we have no additional shells beyond correlated ones,
            # and the data will be copied from corr_shells into shells (see below)
            # number of corr. shells (e.g. Fe d, Ce f) in the unit cell,
            n_corr_shells = int(next(R))
            # now read the information about the correlated shells (atom, sort,
            # l, dim, SO flag, irep):
            corr_shells = [{name: int(val) for name, val in zip(
                corr_shell_entries, R)} for icrsh in range(n_corr_shells)]
            try:
                self.fermi_energy = float(next(R))
            except:
                self.fermi_energy = 0.
        except StopIteration:  # a more explicit error if the file is corrupted.
            mpi.report(self._name + ": reading input file %s failed!" %
                       self.inp_file)
        # close the input file
        R.close()

        # Set or derive some quantities
        # Wannier90 does not use symmetries to reduce the k-points
        # the following might change in future versions
        symm_op = 0
        # copy corr_shells into shells (see above)
        n_shells = n_corr_shells
        shells = []
        for ish in range(n_shells):
            shells.append({key: corr_shells[ish].get(
                key, None) for key in shell_entries})
        ###
        SP = 0                          # NO spin-polarised calculations for now
        SO = 0                          # NO spin-orbit calculation for now
        charge_below = 0                # total charge below energy window NOT used for now
        energy_unit = 1.0               # should be understood as eV units
        ###
        # this is more general
        n_spin = SP + 1 - SO
        dim_corr_shells = sum([sh['dim'] for sh in corr_shells])
        mpi.report(
            "Total number of WFs expected in the correlated shells: %d" % dim_corr_shells)

        # determine the number of inequivalent correlated shells and maps,
        # needed for further processing
        n_inequiv_shells, corr_to_inequiv, inequiv_to_corr = ConverterTools.det_shell_equivalence(
            self, corr_shells)
        mpi.report("Number of inequivalent shells: %d" % n_inequiv_shells)
        mpi.report("Shell representatives: " + format(inequiv_to_corr))
        shells_map = [inequiv_to_corr[corr_to_inequiv[ish]]
                      for ish in range(n_corr_shells)]
        mpi.report("Mapping: " + format(shells_map))
        mpi.report("Subtracting %f eV from the Fermi level." % self.fermi_energy)

        # build the k-point mesh, if its size was given on input (kmesh_mode >= 0),
        # otherwise it is built according to the data in the hr file (see
        # below)
        if kmesh_mode >= 0:
            n_k, k_mesh, bz_weights = self.kmesh_build(nki, kmesh_mode)
            self.n_k = n_k
            self.k_mesh = k_mesh

        # not used in this version: reset to dummy values?
        n_reps = [1 for i in range(n_inequiv_shells)]
        dim_reps = [0 for i in range(n_inequiv_shells)]
        T = []
        for ish in range(n_inequiv_shells):
            ll = 2 * corr_shells[inequiv_to_corr[ish]]['l'] + 1
            lmax = ll * (corr_shells[inequiv_to_corr[ish]]['SO'] + 1)
            T.append(numpy.zeros([lmax, lmax], numpy.complex_))

        spin_w90name = ['_up', '_down']
        hamr_full = []

        # TODO: generalise to SP=1 (only partially done)
        rot_mat_time_inv = [0 for i in range(n_corr_shells)]

        # Second, let's read the file containing the Hamiltonian in WF basis
        # produced by Wannier90
        for isp in range(n_spin):
            # begin loop on isp

            # build filename according to wannier90 conventions
            if SP == 1:
                mpi.report(
                    "Reading information for spin component n. %d" % isp)
                hr_file = self.w90_seed + spin_w90name[isp] + '_hr.dat'
            else:
                hr_file = self.w90_seed + '_hr.dat'
            # now grab the data from the H(R) file
            mpi.report(
                "The Hamiltonian in MLWF basis is extracted from %s ..." % hr_file)
            nr, rvec, rdeg, nw, hamr = self.read_wannier90hr(hr_file)
            # number of R vectors, their indices, their degeneracy, number of
            # WFs, H(R)
            mpi.report("... done: %d R vectors, %d WFs found" % (nr, nw))

            if isp == 0:
                # set or check some quantities that must be the same for both
                # spins
                self.nrpt = nr

                # k-point grid: (if not defined before)
                if kmesh_mode == -1:
                    # the size of the k-point mesh is determined from the
                    # largest R vector
                    nki = [2 * rvec[:, idir].max() + 1 for idir in range(3)]
                    # it will be the same as in the win only when nki is odd, because of the
                    # wannier90 convention: if we have nki k-points along the i-th direction,
                    # then we should get 2*(nki/2)+nki%2 R points along that
                    # direction
                    n_k, k_mesh, bz_weights = self.kmesh_build(nki)
                self.n_k = n_k
                self.k_mesh = k_mesh

                # set the R vectors and their degeneracy
                self.rvec = rvec
                self.rdeg = rdeg

                self.nwfs = nw
                # check that the total number of WFs makes sense
                if self.nwfs < dim_corr_shells:
                    mpi.report(
                        "ERROR: number of WFs in the file smaller than number of correlated orbitals!")
                elif self.nwfs > dim_corr_shells:
                    # NOTE: correlated shells must appear before uncorrelated
                    # ones inside the file
                    mpi.report("Number of WFs larger than correlated orbitals:\n" +
                               "WFs from %d to %d treated as uncorrelated" % (dim_corr_shells + 1, self.nwfs))
                else:
                    mpi.report(
                        "Number of WFs equal to number of correlated orbitals")

                # we assume spin up and spin down always have same total number
                # of WFs
                n_orbitals = numpy.ones(
                    [self.n_k, n_spin], numpy.int) * self.nwfs

            else:
                # consistency check between the _up and _down file contents
                if nr != self.nrpt:
                    mpi.report(
                        "Different number of R vectors for spin-up/spin-down!")
                if nw != self.nwfs:
                    mpi.report(
                        "Different number of WFs for spin-up/spin-down!")

            hamr_full.append(hamr)
            # FIXME: when do we actually need deepcopy()?
            # hamr_full.append(deepcopy(hamr))

            for ir in range(nr):
                # checks if the Hamiltonian is real (it should, if
                # wannierisation worked fine)
                if numpy.abs((hamr[ir].imag.max()).max()) > self._w90zero:
                    mpi.report(
                        "H(R) has large complex components at R %d" % ir)
                # copy the R=0 block corresponding to the correlated shells
                # into another variable (needed later for finding rot_mat)
                if rvec[ir, 0] == 0 and rvec[ir, 1] == 0 and rvec[ir, 2] == 0:
                    ham_corr0 = hamr[ir][0:dim_corr_shells, 0:dim_corr_shells]

            # checks if ham0 is Hermitian
            if not numpy.allclose(ham_corr0.transpose().conjugate(), ham_corr0, atol=self._w90zero, rtol=1.e-9):
                raise ValueError("H(R=0) matrix is not Hermitian!")

            # find rot_mat symmetries by diagonalising the on-site Hamiltonian
            # of the first spin
            if isp == 0:
                use_rotations, rot_mat = self.find_rot_mat(
                    n_corr_shells, corr_shells, shells_map, ham_corr0)
            else:
                # consistency check
                use_rotations_, rot_mat_ = self.find_rot_mat(
                    n_corr_shells, corr_shells, shells_map, ham_corr0)
                if (use_rotations and not use_rotations_):
                    mpi.report(
                        "Rotations cannot be used for spin component n. %d" % isp)
                for icrsh in range(n_corr_shells):
                    if not numpy.allclose(rot_mat_[icrsh], rot_mat[icrsh], atol=self._w90zero, rtol=1.e-15):
                        mpi.report(
                            "Rotations for spin component n. %d do not match!" % isp)
        # end loop on isp

        mpi.report("The k-point grid has dimensions: %d, %d, %d" % tuple(nki))
        # if calculations are spin-polarized, then renormalize k-point weights
        if SP == 1:
            bz_weights = 0.5 * bz_weights

        # Third, compute the hoppings in reciprocal space
        hopping = numpy.zeros([self.n_k, n_spin, numpy.max(
            n_orbitals), numpy.max(n_orbitals)], numpy.complex_)
        for isp in range(n_spin):
            # make Fourier transform H(R) -> H(k) : it can be done one spin at
            # a time
            hamk = self.fourier_ham(self.nwfs, hamr_full[isp])
            # copy the H(k) in the right place of hoppings... is there a better
            # way to do this??
            for ik in range(self.n_k):
                #hopping[ik,isp,:,:] = deepcopy(hamk[ik][:,:])*energy_unit
                hopping[ik, isp, :, :] = hamk[ik][:, :] * energy_unit

        # Then, initialise the projectors
        k_dep_projection = 0   # we always have the same number of WFs at each k-point
        proj_mat = numpy.zeros([self.n_k, n_spin, n_corr_shells, max(
            [crsh['dim'] for crsh in corr_shells]), numpy.max(n_orbitals)], numpy.complex_)
        iorb = 0
        # Projectors simply consist in identity matrix blocks selecting those MLWFs that
        # correspond to the specific correlated shell indexed by icrsh.
        # NOTE: we assume that the correlated orbitals appear at the beginning of the H(R)
        # file and that the ordering of MLWFs matches the corr_shell info from
        # the input.
        for icrsh in range(n_corr_shells):
            norb = corr_shells[icrsh]['dim']
            proj_mat[:, :, icrsh, 0:norb, iorb:iorb +
                     norb] = numpy.identity(norb, numpy.complex_)
            iorb += norb

        # Finally, save all required data into the HDF archive:
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.dft_subgrp in ar):
                ar.create_group(self.dft_subgrp)
            # The subgroup containing the data. If it does not exist, it is
            # created. If it exists, the data is overwritten!
            things_to_save = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                          'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                          'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights', 'hopping',
                          'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
            for it in things_to_save:
                ar[self.dft_subgrp][it] = locals()[it]

    def read_wannier90hr(self, hr_filename="wannier_hr.dat"):
        """
        Method for reading the seedname_hr.dat file produced by Wannier90 (http://wannier.org)

        Parameters
        ----------
        hr_filename : string
            full name of the H(R) file produced by Wannier90 (usually seedname_hr.dat)

        Returns
        -------
        nrpt : integer
            number of R vectors found in the file
        rvec_idx : numpy.array of integers
            Miller indices of the R vectors
        rvec_deg : numpy.array of floats
            weight of the R vectors
        num_wf : integer
            number of Wannier functions found
        h_of_r : list of numpy.array
            <w_i|H(R)|w_j> = Hamilonian matrix elements in the Wannier basis

        """

        # Read only from the master node
        if not (mpi.is_master_node()):
            return

        try:
            with open(hr_filename, "r") as hr_filedesc:
                hr_data = hr_filedesc.readlines()
                hr_filedesc.close()
        except IOError:
            mpi.report("The file %s could not be read!" % hr_filename)

        mpi.report("Reading %s..." % hr_filename + hr_data[0])

        try:
            # reads number of Wannier functions per spin
            num_wf = int(hr_data[1])
            nrpt = int(hr_data[2])
        except ValueError:
            mpi.report("Could not read number of WFs or R vectors")

        # allocate arrays to save the R vector indexes and degeneracies and the
        # Hamiltonian
        rvec_idx = numpy.zeros((nrpt, 3), dtype=int)
        rvec_deg = numpy.zeros(nrpt, dtype=int)
        h_of_r = [numpy.zeros((num_wf, num_wf), dtype=numpy.complex_)
                  for n in range(nrpt)]

        # variable currpos points to the current line in the file
        currpos = 2
        try:
            ir = 0
            # read the degeneracy of the R vectors (needed for the Fourier
            # transform)
            while ir < nrpt:
                currpos += 1
                for x in hr_data[currpos].split():
                    if ir >= nrpt:
                        raise IndexError("wrong number of R vectors??")
                    rvec_deg[ir] = int(x)
                    ir += 1
            # for each direct lattice vector R read the block of the
            # Hamiltonian H(R)
            for ir, jj, ii in product(list(range(nrpt)), list(range(num_wf)), list(range(num_wf))):
                # advance one line, split the line into tokens
                currpos += 1
                cline = hr_data[currpos].split()
                # check if the orbital indexes in the file make sense
                if int(cline[3]) != ii + 1 or int(cline[4]) != jj + 1:
                    mpi.report(
                        "Inconsistent indices at %s%s of R n. %s" % (ii, jj, ir))
                rcurr = numpy.array(
                    [int(cline[0]), int(cline[1]), int(cline[2])])
                if ii == 0 and jj == 0:
                    rvec_idx[ir] = rcurr
                    rprec = rcurr
                else:
                    # check if the vector indices are consistent
                    if not numpy.array_equal(rcurr, rprec):
                        mpi.report(
                            "Inconsistent indices for R vector n. %s" % ir)

                # fill h_of_r with the matrix elements of the Hamiltonian
                if not numpy.any(rcurr) and ii == jj:
                    h_of_r[ir][ii, jj] = complex(float(cline[5]) - self.fermi_energy, float(cline[6]))
                else:
                    h_of_r[ir][ii, jj] = complex(float(cline[5]), float(cline[6]))

        except ValueError:
            mpi.report("Wrong data or structure in file %s" % hr_filename)

        # return the data into variables
        return nrpt, rvec_idx, rvec_deg, num_wf, h_of_r

    def find_rot_mat(self, n_sh, sh_lst, sh_map, ham0):
        """
        Method for finding the matrices that bring from local to global coordinate systems
        (and viceversa), based on the eigenvalues of H(R=0)

        Parameters
        ----------
        n_sh : integer
            number of shells
        sh_lst : list of shells-type dictionaries
            contains the shells (could be correlated or not)
        sh_map : list of integers
            mapping between shells
        ham0 : numpy.array of floats
            local Hamiltonian matrix elements

        Returns
        -------
        istatus : integer
            if 0, something failed in the construction of the matrices
        rot_mat : list of numpy.array
            rotation matrix for each of the shell

        """

        # initialize the rotation matrices to identities
        rot_mat = [numpy.identity(sh_lst[ish]['dim'], dtype=complex)
                   for ish in range(n_sh)]
        istatus = 0

        hs = ham0.shape
        if hs[0] != hs[1] or hs[0] != sum([sh['dim'] for sh in sh_lst]):
            mpi.report(
                "find_rot_mat: wrong block structure of input Hamiltonian!")
            istatus = 0
            # this error will lead into troubles later... early return
            return istatus, rot_mat

        # TODO: better handling of degenerate eigenvalue case
        eigval_lst = []
        eigvec_lst = []
        iwf = 0
        # loop over shells
        for ish in range(n_sh):
            # nw = number of orbitals in this shell
            nw = sh_lst[ish]["dim"]
            # diagonalize the sub-block of H(0) corresponding to this shell
            eigval, eigvec = numpy.linalg.eigh(
                ham0[iwf:iwf + nw, iwf:iwf + nw])
            # find the indices sorting the eigenvalues in ascending order
            eigsrt = eigval[0:nw].argsort()
            # order eigenvalues and eigenvectors and save in a list
            eigval_lst.append(eigval[eigsrt])
            eigvec_lst.append(eigvec[eigsrt])
            iwf += nw
            # TODO: better handling of degenerate eigenvalue case
            if sh_map[ish] != ish:  # issue warning only when there are equivalent shells
                for i in range(nw):
                    for j in range(i + 1, nw):
                        if (abs(eigval[j] - eigval[i]) < self._w90zero):
                            mpi.report("WARNING: degenerate eigenvalue of H(0) detected for shell %d: " % (ish) +
                                       "global-to-local transformation might not work!")

        for ish in range(n_sh):
            try:
                # build rotation matrices by combining the unitary
                # transformations that diagonalize H(0)
                rot_mat[ish] = numpy.dot(eigvec_lst[ish], eigvec_lst[
                                         sh_map[ish]].conjugate().transpose())
            except ValueError:
                mpi.report(
                    "Global-to-local rotation matrices cannot be constructed!")

            istatus = 1
            # check that eigenvalues are the same (within accuracy) for
            # equivalent shells
            if not numpy.allclose(eigval_lst[ish], eigval_lst[sh_map[ish]], atol=self._w90zero, rtol=1.e-15):
                mpi.report(
                    "ERROR: eigenvalue mismatch between equivalent shells! %d" % ish)
                eigval_diff = eigval_lst[ish] - eigval_lst[sh_map[ish]]
                mpi.report("Eigenvalue difference: " + format(eigval_diff))
                istatus = 0

        # TODO: add additional consistency check on rot_mat matrices?

        return istatus, rot_mat

    def kmesh_build(self, msize=None, mmode=0):
        """
        Method for the generation of the k-point mesh.
        Right now it only supports the option for generating a full grid containing k=0,0,0.

        Parameters
        ----------
        msize : list of 3 integers
            the dimensions of the mesh
        mmode : integer
            mesh generation mode (right now, only full grid available)

        Returns
        -------
        nkpt : integer
            total number of k-points in the mesh
        k_mesh : numpy.array[nkpt,3] of floats
            the coordinates of all k-points
        wk : numpy.array[nkpt] of floats
            the weight of each k-point

        """

        if mmode != 0:
            raise ValueError("Mesh generation mode not supported: %s" % mmode)

        # a regular mesh including Gamma point
        # total number of k-points
        nkpt = msize[0] * msize[1] * msize[2]
        kmesh = numpy.zeros((nkpt, 3), dtype=float)
        ii = 0
        for ix, iy, iz in product(list(range(msize[0])), list(range(msize[1])), list(range(msize[2]))):
            kmesh[ii, :] = [float(ix) / msize[0], float(iy) /
                            msize[1], float(iz) / msize[2]]
            ii += 1
        # weight is equal for all k-points because wannier90 uses uniform grid on whole BZ
        # (normalization is always 1 and takes into account spin degeneracy)
        wk = numpy.ones([nkpt], dtype=float) / float(nkpt)

        return nkpt, kmesh, wk

    def fourier_ham(self, norb, h_of_r):
        """
        Method for obtaining H(k) from H(R) via Fourier transform
        The R vectors and k-point mesh are read from global module variables

        Parameters
        ----------
        norb : integer
            number of orbitals
        h_of_r : list of numpy.array[norb,norb]
            Hamiltonian H(R) in Wannier basis

        Returns
        -------
        h_of_k : list of numpy.array[norb,norb]
            transformed Hamiltonian H(k) in Wannier basis

        """

        twopi = 2 * numpy.pi
        h_of_k = [numpy.zeros((norb, norb), dtype=numpy.complex_)
                  for ik in range(self.n_k)]
        ridx = numpy.array(list(range(self.nrpt)))
        for ik, ir in product(list(range(self.n_k)), ridx):
            rdotk = twopi * numpy.dot(self.k_mesh[ik], self.rvec[ir])
            factor = (math.cos(rdotk) + 1j * math.sin(rdotk)) / \
                float(self.rdeg[ir])
            h_of_k[ik][:, :] += factor * h_of_r[ir][:, :]

        return h_of_k
