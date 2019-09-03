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

import os

from .base import LatticeModel
from .tools import set_nk
from pytriqs.archive.hdf_archive import HDFArchive

class ExternalModel(LatticeModel):
    """
    Prepare DFT data externally. This class only checks an existing file and data in it.
    """

    def __init__(self, params):
        super(ExternalModel, self).__init__(params)

        from pytriqs.archive import HDFArchive

        seedname = self._params["model"]["seedname"]
        h5_file = seedname+'.h5'

        # check if h5 file already exists
        try:
            assert os.path.exists(h5_file)
            with HDFArchive(h5_file, 'r') as ar:
                assert 'dft_input' in ar
        except:
            raise Exception("Prepare, in advance, '%s' file which stores DFT data in 'dft_input' subgroup" % h5_file)

        # set nkdiv
        self._nkdiv = set_nk(params["model"]["nk"],
                             params["model"]["nk0"],
                             params["model"]["nk1"],
                             params["model"]["nk2"])

        # Set [model][norb_inequiv_sh], which is necessary for generate_umat
        with HDFArchive(h5_file, 'r') as f:
            n_inequiv_shells = f["dft_input"]["n_inequiv_shells"]
            corr_shells = f["dft_input"]["corr_shells"]
            inequiv_to_corr = f["dft_input"]["inequiv_to_corr"]
            norb_inequiv_sh = [corr_shells[icsh]['dim'] for icsh in inequiv_to_corr]
            try:
                assert len(norb_inequiv_sh) == n_inequiv_shells
            except:
                print("ExternalModel.__ini__: failed in setting 'norb_inequiv_sh'")
                print(" n_inequiv_shells =", n_inequiv_shells)
                print(" inequiv_to_corr =", inequiv_to_corr)
                print(" corr_shells =", corr_shells)
                print(" norb_inequiv_sh =", norb_inequiv_sh)
                exit(1)
            params['model']['norb_inequiv_sh'] = norb_inequiv_sh

    @classmethod
    def name(self):
        return 'external'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        pass
