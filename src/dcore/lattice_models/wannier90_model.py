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
from h5 import HDFArchive

from .base import LatticeModel
from .tools import set_nk
from dcore.converters.wannier90 import Wannier90Converter

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
    corr_to_inequiv = params['model']['corr_to_inequiv']
    if params['model']['spin_orbit']:
        dim_Hk = 2 * params['model']['norb_corr_sh']
    else:
        dim_Hk = params['model']['norb_corr_sh']

    nk0, nk1, nk2 = nkdiv

    print("\n    nk0 = {0}".format(nk0))
    print("    nk1 = {0}".format(nk1))
    print("    nk2 = {0}".format(nk2))
    print("    ncor = {0}".format(ncor))
    print("    n_inequiv_sh = {0}".format(params['model']['n_inequiv_shells']))
    for i in range(ncor):
        assert corr_to_inequiv[i] >= 0
        print("    dim[{0}], equiv[{0}] = {1}, {2}".format(i, dim_Hk[i], corr_to_inequiv[i]))
    print("")

    #
    # Generate file for Wannier90-Converter
    #
    print("0 {0} {1} {2}".format(nk0, nk1, nk2), file=f)
    print(params["model"]["nelec"], file=f)
    print(ncor, file=f)
    for i in range(ncor):
        print("{0} {1} {2} {3} 0 0".format(i, corr_to_inequiv[i], 0, dim_Hk[i]), file=f)


class Wannier90Model(LatticeModel):
    def __init__(self, params):
        super(Wannier90Model, self).__init__(params)

        self._nkdiv = set_nk(params["model"]["nk"],
                             params["model"]["nk0"],
                             params["model"]["nk1"],
                             params["model"]["nk2"])
        self._spin_orbit = params['model']['spin_orbit']


    @classmethod
    def name(self):
        return 'wannier90'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        #
        #
        p = self._params
        seedname = self._params["model"]["seedname"]
        with open(seedname+'.inp', 'w') as f:
            _generate_w90_converter_input(self.nkdiv(), p, f)

        # Convert General-Hk to SumDFT-HDF5 format
        converter = Wannier90Converter(seedname=seedname)
        converter.convert_dft_input()

        if p["model"]["spin_orbit"]:
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
                # (n_k, nb, n_corr, nso, orb, spin) => (n_k, nb, n_corr, nso, spin, orb)
                assert max_dim_sh//2 > 0
                proj_mat = proj_mat.reshape((n_k, nb, n_corr, max_dim_sh, max_n_orb//2, 2))
                proj_mat = proj_mat.transpose((0, 1, 2, 3, 5, 4))
                proj_mat = proj_mat.reshape((n_k, 1, n_corr, max_dim_sh, max_n_orb))
                f['dft_input']['proj_mat'] = proj_mat


    def write_dft_band_input_data(self, params, kvec, bands_data='dft_bands_input'):
        """

        Returns
        -------
        hopping : complex
            k-dependent one-body Hamiltonian
        n_orbitals : integer
            Number of orbitals at each k. It does not depend on k
        proj_mat : complex
            Projection onto each correlated orbitals
        band_data : str
            Where the data is to be written
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
            offset += dim_Hk[icor]

        # Make proj_mat compatible with DCore's block structure
        if spin_orbit:
            proj_mat = proj_mat.reshape((n_k, nblock, ncor, numpy.max(dim_Hk)//2, 2, nwan//2, 2))
            proj_mat = proj_mat.transpose((0, 1, 2, 4, 3, 6, 5))
            proj_mat = proj_mat.reshape((n_k, nblock, ncor, numpy.max(dim_Hk), nwan))

        #
        # Output them into seedname.h5
        #
        with HDFArchive(seedname + '.h5', 'a') as f:
            if not (bands_data in f):
                f.create_group(bands_data)
            f[bands_data]['hopping'] = hopping
            f[bands_data]['n_k'] = n_k
            f[bands_data]['n_orbitals'] = numpy.full((n_k, nblock), nwan, dtype=int)
            f[bands_data]['proj_mat'] = proj_mat
        print('    Done')
