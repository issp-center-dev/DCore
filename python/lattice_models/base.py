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

    @classmethod
    def is_Hk_supported(cls):
        return True

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

    def write_dft_band_input_data(self, params, kvec, bands_data='dft_bands_input'):
        """

        Parameters
        ----------
        params
        kvec: float of (*, 3)
            a list of k points. Each element should be between 0 and 2*pi.
        bands_data

        Returns
        -------

        """
        pass

