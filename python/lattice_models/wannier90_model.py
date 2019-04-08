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

# DO NOT IMPORT GLOBALLY ANY MODULE DEPENDING ON MPI
#import os
import re
#import sys
#import numpy
#import scipy
#from itertools import product
from pytriqs.archive.hdf_archive import HDFArchive

from .base import LatticeModel

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
    ncor = params["model"]["ncor"]
    #
    norb_list = re.findall(r'\d+', params["model"]["norb"])
    norb = [int(norb_list[icor]) for icor in range(ncor)]
    #
    if params["model"]["equiv"] == "None":
        equiv = [icor for icor in range(ncor)]
    else:
        equiv_list = re.findall(r'[^\s,]+', params["model"]["equiv"])
        equiv_str_list = []
        equiv_index = 0
        equiv = [0] * ncor
        for icor in range(ncor):
            if equiv_list[icor] in equiv_str_list:
                # Defined before
                equiv[icor] = equiv_str_list.index(equiv_list[icor])
            else:
                # New one
                equiv_str_list.append(equiv_list[icor])
                equiv[icor] = equiv_index
                equiv_index += 1
    nk0, nk1, nk2 = nkdiv

    print("\n    nk0 = {0}".format(nk0))
    print("    nk1 = {0}".format(nk1))
    print("    nk2 = {0}".format(nk2))
    print("    ncor = {0}".format(ncor))
    for i in range(ncor):
        assert equiv[i] >= 0
        print("    norb[{0}], equiv[{0}] = {1}, {2}".format(i, norb[i], equiv[i]))
    print("")
    #
    # Generate file for Wannier90-Converter
    #
    print("0 {0} {1} {2}".format(nk0, nk1, nk2), file=f)
    print(params["model"]["nelec"], file=f)
    print(ncor, file=f)
    for i in range(ncor):
        print("{0} {1} {2} {3} 0 0".format(i, equiv[i], 0, norb[i]), file=f)


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

    @classmethod
    def name(self):
        return 'wannier90'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        from ..converters.wannier90_converter import Wannier90Converter

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

        if p["model"]["spin_orbit"] or p["model"]["non_colinear"]:
            if p['model']['non_colinear']:
                from ..manip_database import turn_on_spin_orbit
                print('')
                print('Turning on spin_orbit...')
                turn_on_spin_orbit(seedname + '.h5', seedname + '.h5')
            else:
                with HDFArchive(seedname + '.h5', 'a') as f:
                    f["dft_input"]["SP"] = 1
                    f["dft_input"]["SO"] = 1

                    corr_shells = f["dft_input"]["corr_shells"]
                    for icor in range(p["model"]['ncor']):
                        corr_shells[icor]["SO"] = 1
                    f["dft_input"]["corr_shells"] = corr_shells

