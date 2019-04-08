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

from pytriqs.archive.hdf_archive import HDFArchive


class LatticeModel(object):
    """
    Base class for H(k) model
    """
    def __init__(self, params):
        self._nkdiv = (1, 1, 1)
        self._params = params

    @classmethod
    def name(cls):
        return 'lattice_base_model'

    def nkdiv(self):
        return self._nkdiv

    def Hk(self, kvec):
        """

        kvec is a 1D array-like object of length 3 (with 2*pi)

        :return:
            spin_orbit == True:
                Hk for ud sector: complex(2*num_orb, 2*num_orb)
            spin_orbit == False:
                Hk for up and down sectors: (complex(num_orb, num_orb), complex(num_orb, num_orb))
        """
        pass

    def generate_model_file(self):
        pass

    def write_dft_band_input_data(self, params, kvec):
        pass


def write_dft_bands_input_data(seedname, params, n_k, kvec, lattice_model):
    """
    Write DFT band input data into a HDF5 file
    :param seedname:
    :param params: runtime parameters
    :param n_k:  Number of k points
    :param kvec:  2D array
    :param lattice_model: object of LatticeModel
    :return:  None
    """
    assert kvec.shape[0] == n_k

    hopping, n_orbitals, proj_mat = lattice_model.generate_dft_band_input_data(params, kvec)

    #
    # Output them into seedname.h5
    #
    with HDFArchive(seedname+'.h5', 'a') as f:
        if not ('dft_bands_input' in f):
            f.create_group('dft_bands_input')
        f['dft_bands_input']['hopping'] = hopping
        f['dft_bands_input']['n_k'] = n_k
        f['dft_bands_input']['n_orbitals'] = n_orbitals
        f['dft_bands_input']['proj_mat'] = proj_mat
    print('    Done')

