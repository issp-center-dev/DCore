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


import os
import numpy
from warnings import warn

from .base import LatticeModel
from .tools import set_nk, XNode
from h5 import HDFArchive


class ExternalModel(LatticeModel):
    """
    Prepare DFT data externally. This class only checks an existing file and data in it.
    """

    def __init__(self, params):
        super(ExternalModel, self).__init__(params)

        self._seedname = self._params["model"]["seedname"]
        h5_file = self._seedname+'.h5'

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
                print("ExternalModel.__init__: failed in setting 'norb_inequiv_sh'")
                print(" n_inequiv_shells =", n_inequiv_shells)
                print(" inequiv_to_corr =", inequiv_to_corr)
                print(" corr_shells =", corr_shells)
                print(" norb_inequiv_sh =", norb_inequiv_sh)
                exit(1)
            params['model']['norb_inequiv_sh'] = norb_inequiv_sh

    @classmethod
    def name(cls):
        return 'external'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        pass

    def is_Hk_supported(cls):
        return False

    def generate_Hk_path(self, params):
        """
        Override the original method. It is assumed that H(k) data are already stored in seedname.h5.
        """

        # check h5 file
        with HDFArchive(self._seedname + '.h5', 'r') as f:
            if not 'dft_bands_input' in f:
                warn("  seedname.h5/dft_bands_input should be prepared in advance in lattice='external'")
                return None, None

            n_k = f["dft_bands_input"]["n_k"]

        # generate x values
        xk = []
        xnode = []
        file_kpath = "kpath.in"
        if not os.path.exists(file_kpath):
            warn("  '%s' not found" %file_kpath)
            return None, None
        with open(file_kpath, "r") as f:
            # -40 -40 40 80 2.31040 Z
            for line in f:
                array = line.split()
                assert len(array) >= 5
                x = float(array[4])
                xk.append(x)
                if len(array) >= 6:
                    xnode.append(XNode(x = x, label = array[5]))

        assert len(xk) == n_k
        return numpy.array(xk), xnode
