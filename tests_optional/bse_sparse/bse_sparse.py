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

from dcore.dcore_pre import dcore_pre
from dcore.dcore import dcore
from dcore.dcore_bse import dcore_bse
from bse_tools.h5bse import h5BSE
import shutil
import numpy
import sys
import unittest




class TestMethods(unittest.TestCase):

    def setUp(self):
        # reference data (non-sparse sampling)
        bse = h5BSE("dmft_bse.h5.ref")
        self.xloc_ref = bse.get(('X_loc', 0))

    def _fit_and_compare(self, postfix):

        dcore_bse("sparse_fit.ini" + postfix)
        sys.stdout.flush()

        # interpolated data
        bse = h5BSE("dmft_bse.h5" + postfix)
        xloc_fit = bse.get(('X_loc', 0))

        # compare
        for key in xloc_fit:
            # for test
            print("distance:", numpy.linalg.norm(xloc_fit[key] - self.xloc_ref[key]))
            sys.stdout.flush()

            # compare
            numpy.testing.assert_allclose(xloc_fit[key], self.xloc_ref[key], atol=1e-2)
            # self.assertTrue(numpy.allclose(xloc_fit[key], self.xloc_ref[key], rtol=1e-2))
            sys.stdout.flush()

    def test_fit_Lambda100(self):
        self._fit_and_compare(".Lambda100")


if __name__ == '__main__':
    input_ini = "sparse_fit.ini.ref"
    dcore_pre(input_ini)
    dcore(input_ini)
    dcore_bse(input_ini)

    unittest.main()
