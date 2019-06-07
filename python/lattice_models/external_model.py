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

    @classmethod
    def name(self):
        return 'external'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        pass
