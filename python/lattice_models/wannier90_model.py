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

import re
import numpy
from pytriqs.archive.hdf_archive import HDFArchive

from .base import LatticeModel
from ..converters.wannier90_converter import Wannier90Converter

def _generate_w90_converter_input(nkdiv, params, f):
    """
    Compute hopping etc. for A(k,w) of Wannier90

    Parameters
    ----------
    nkdiv : (int, int, int)
        Number of k divisions
    params : dictionary
        Input parameters
    f : File stream
        file
    """

    #
    # Number of orbitals in each shell
    # Equivalence of each shell
    #
    ncor = params['model']['ncor']
    equiv_sh = params['model']['equiv_sh']
    if params['model']['spin_orbit']:
        dim_Hk = 2 * params['model']['norb_corr_sh']
    else:
        dim_Hk = params['model']['norb_corr_sh']

    #
    #norb_list = re.findall(r'\d+', params["model"]["norb"])
    #norb = [int(norb_list[icor]) for icor in range(ncor)]
    #
    nk0, nk1, nk2 = nkdiv

    print("\n    nk0 = {0}".format(nk0))
    print("    nk1 = {0}".format(nk1))
    print("    nk2 = {0}".format(nk2))
    print("    ncor = {0}".format(ncor))
    for i in range(ncor):
        assert equiv_sh[i] >= 0
        print("    dim[{0}], equiv[{0}] = {1}, {2}".format(i, dim_Hk[i], equiv_sh[i]))
    print("")

    #
    # Generate file for Wannier90-Converter
    #
    print("0 {0} {1} {2}".format(nk0, nk1, nk2), file=f)
    print(params["model"]["nelec"], file=f)
    print(ncor, file=f)
    for i in range(ncor):
        print("{0} {1} {2} {3} 0 0".format(i, equiv_sh[i], 0, dim_Hk[i]), file=f)


def _set_nk(nk, nk0, nk1, nk2):
    if abs(nk0) + abs(nk1) + abs(nk2) == 0:
        # If one of nk0, nk1 and nk2 are set, use nk.
        nk0 = nk
        nk1 = nk
        nk2 = nk
    elif abs(nk0) + abs(nk1) + abs(nk2) > 0:
        # if any of nk0, nk1 and nk2 are set, use them.
        if nk0 * nk1 * nk2 == 0:
            raise RuntimeError("Some of nk0, nk1 and nk2 are zero!")
    return nk0, nk1, nk2

class Wannier90Model(LatticeModel):
    def __init__(self, params):
        super(Wannier90Model, self).__init__(params)

        self._nkdiv = _set_nk(params["system"]["nk"],
                             params["system"]["nk0"],
                             params["system"]["nk1"],
                             params["system"]["nk2"])
        self._spin_orbit = params['model']['spin_orbit']

        #self._norb_list = params['model']['norb_list']
        #norb = [int(norb_list[icor]) for icor in range(ncor)]

    @classmethod
    def name(self):
        return 'wannier90'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        #
        # non_colinear flag is used only for the case that COLINEAR DFT calculation
        #
        p = self._params
        seedname = self._params["model"]["seedname"]
        if p["model"]["spin_orbit"] and p["model"]["non_colinear"]:
            raise RuntimeError("The options spin_orbit and non_colinear are not compatible!")
        with open(seedname+'.inp', 'w') as f:
            _generate_w90_converter_input(self.nkdiv(), p, f)

        # Convert General-Hk to SumDFT-HDF5 format
        converter = Wannier90Converter(seedname=seedname)
        converter.convert_dft_input()

        if p['model']['non_colinear']:
            from ..manip_database import turn_on_spin_orbit
            print('')
            print('Turning on spin_orbit...')
            turn_on_spin_orbit(seedname + '.h5', seedname + '.h5')
        elif p["model"]["spin_orbit"]:
            with HDFArchive(seedname + '.h5', 'a') as f:
                f["dft_input"]["SP"] = 1
                f["dft_input"]["SO"] = 1

                corr_shells = f["dft_input"]["corr_shells"]
                for icor in range(p["model"]['ncor']):
                    corr_shells[icor]["SO"] = 1
                f["dft_input"]["corr_shells"] = corr_shells

                # Make projectors compatible with DCore's block structure
                # proj_mat (n_k, nb, n_corr, max_dim_sh, max_n_orb)
                proj_mat = f['dft_input']['proj_mat']
                n_k, nb, n_corr, max_dim_sh, max_n_orb = proj_mat.shape
                assert nb == 1
                # (n_k, nb, n_corr, orb, spin, orb, spin) => (n_k, nb, n_corr, spin, orb, spin, orb)
                proj_mat = proj_mat.reshape((n_k, nb, n_corr, max_dim_sh//2, 2, max_n_orb//2, 2))
                proj_mat = proj_mat.transpose((0, 1, 2, 4, 3, 6, 5))
                f['dft_input']['proj_mat'] = proj_mat.reshape((n_k, 1, n_corr, max_dim_sh, max_n_orb))


    def write_dft_band_input_data(self, params, kvec):
        """

        Returns
        -------
        hopping : complex
            k-dependent one-body Hamiltonian
        n_orbitals : integer
            Number of orbitals at each k. It does not depend on k
        proj_mat : complex
            Projection onto each correlated orbitals
        """
        n_k = kvec.shape[0]
        assert kvec.shape[1] == 3

        spin_orbit = params['model']['spin_orbit']

        norb_sh = numpy.array(params['model']['norb_corr_sh'])
        n_spin_orb_sh = 2 * norb_sh
        ncor = params['model']['ncor']
        seedname = params['model']['seedname']

        # Print 2 * norb for each shell
        print("               ncor = ", ncor)
        for i in range(ncor):
            print("     num_spin_orb[{0}] = {1}".format(i, n_spin_orb_sh[i]))

        #
        # Read hopping in the real space from the Wannier90 output
        #  FIXME: Wannier90Converter may be replaced by an original implementation
        #
        # if spin_orbit == True, nwan must be twice of the number of orbitals in the correlated shells.
        # Otherwise, nwan must be the number of orbitals in the correlated shells.
        #
        w90c = Wannier90Converter(seedname=seedname)
        nr, rvec, rdeg, nwan, hamr = w90c.read_wannier90hr(seedname + "_hr.dat")

        # Number of physical orbitals in the model (including non-interacting ones)
        #if spin_orbit:
            #norb_tot = nwan//2
        #else:
            #norb_tot = nwan

        #
        # Fourier transformation of the one-body Hamiltonian
        #
        nblock = 1
        hopping = numpy.zeros((n_k, nblock, nwan, nwan), complex)
        for ik in range(n_k):
            for ir in range(nr):
                rdotk = numpy.dot(kvec[ik, :], rvec[ir, :])
                factor = (numpy.cos(rdotk) + 1j * numpy.sin(rdotk)) / float(rdeg[ir])
                hopping[ik, 0, :, :] += factor * hamr[ir][:, :]

        #
        # proj_mat is (norb*norb) identities at each correlation shell
        #
        offset = 0
        dim_Hk = n_spin_orb_sh if spin_orbit else norb_sh
        proj_mat = numpy.zeros([n_k, nblock, ncor, numpy.max(dim_Hk), nwan], complex)
        for icor in range(ncor):
            proj_mat[:, :, icor, 0:dim_Hk[icor], offset:offset+dim_Hk[icor]] = numpy.identity(dim_Hk[icor])
            offset += n_spin_orb_sh[icor]
        # Make proj_mat compatible with DCore's block structure
        proj_mat = proj_mat.reshape((n_k, nblock, ncor, numpy.max(dim_Hk)//2, 2, nwan//2, 2))
        proj_mat = proj_mat.transpose((0, 1, 2, 4, 3, 6, 5))
        proj_mat = proj_mat.reshape((n_k, nblock, ncor, numpy.max(dim_Hk), nwan))

        #
        # Output them into seedname.h5
        #
        with HDFArchive(seedname + '.h5', 'a') as f:
            if not ('dft_bands_input' in f):
                f.create_group('dft_bands_input')
            f['dft_bands_input']['hopping'] = hopping
            f['dft_bands_input']['n_k'] = n_k
            f['dft_bands_input']['n_orbitals'] = numpy.full((n_k, nblock), nwan, dtype=int)
            f['dft_bands_input']['proj_mat'] = proj_mat
        print('    Done')

        #
        # Extend spin block if needed
        #
        if params['model']['non_colinear']:
            from ..manip_database import turn_on_spin_orbit
            print('')
            print('Turning on non_collinear...')
            turn_on_spin_orbit(seedname + '.h5', seedname + '.h5', update_dft_input=False, update_dft_bands_input=True)
