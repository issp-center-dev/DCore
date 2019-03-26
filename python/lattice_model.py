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
import os
import re
import numpy
from warnings import warn


def create_lattice_model(params):
    model_name =  params["model"]["lattice"]

    for c in all_lattice_models:
        if model_name == c.name():
            return c(params)

    raise RuntimeError('Unknown lattice model name: ' + model_name)


def _generate_nnn_lattice_model(spatial_dim, nelec, norb, t, tp, nk, f):
    """
    Compute a list of H(k) for lattice model with t and t'

    Parameters
    ----------
    spatial_dim : int
            Spatial dimension
    nelec : int
            Number of electrons
    norb : int
            Number of orbitals
    t : float
            Nearest neighbor hopping
    tp : float
            Next nearest neighbor hopping
    nk : int
            Number of k divisions along each axis
    f : file
            Output file
    """

    assert spatial_dim > 0 and spatial_dim <= 3

    print("nk = ", nk, spatial_dim)

    nkbz = nk**spatial_dim

    print("\n    Total number of k =", str(nkbz))

    #
    # Write General-Hk formatted file
    #
    print(nkbz, file=f)
    print(nelec, file=f)
    print("1", file=f)
    print("0 0 {0} {1}".format(0, norb), file=f)
    print("1", file=f)
    print("0 0 {0} {1} 0 0".format(0, norb), file=f)
    print("1 {0}".format(norb), file=f)

    #
    # Energy band
    #
    kvec = [0.0, 0.0, 0.0]
    if spatial_dim == 1:
        for i0 in range(nk):
            kvec[0] = 2.0 * numpy.pi * float(i0) / float(nk)
            ek = 2.0*t*numpy.cos(kvec[0]) + 2*tp*numpy.cos(2.0*kvec[0])
            for iorb in range(norb):
                for jorb in range(norb):
                    if iorb == jorb:
                        print("{0}".format(ek), file=f)  # Real part
                    else:
                        print("0.0", file=f)  # Real part
            for iorb in range(norb*norb):
                print("0.0", file=f)  # Imaginary part
    elif spatial_dim == 2:
        for i0 in range(nk):
            kvec[0] = 2.0 * numpy.pi * float(i0) / float(nk)
            for i1 in range(nk):
                kvec[1] = 2.0 * numpy.pi * float(i1) / float(nk)
                ek = 2.0*t*(numpy.cos(kvec[0]) + numpy.cos(kvec[1])) \
                     + 2.0*tp*(numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1]))
                for iorb in range(norb):
                    for jorb in range(norb):
                        if iorb == jorb:
                            print("{0}".format(ek), file=f)  # Real part
                        else:
                            print("0.0", file=f)  # Real part
                for iorb in range(norb*norb):
                    print("0.0", file=f)  # Imaginary part
    elif spatial_dim == 3:
        for i0 in range(nk):
            kvec[0] = 2.0 * numpy.pi * float(i0) / float(nk)
            for i1 in range(nk):
                kvec[1] = 2.0 * numpy.pi * float(i1) / float(nk)
                for i2 in range(nk):
                    kvec[2] = 2.0 * numpy.pi * float(i2) / float(nk)
                    ek = 2*t*(numpy.cos(kvec[0]) + numpy.cos(kvec[1]) + numpy.cos(kvec[2])) \
                         + 2*tp*(numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1])
                                 + numpy.cos(kvec[1] + kvec[2]) + numpy.cos(kvec[1] - kvec[2])
                                 + numpy.cos(kvec[2] + kvec[0]) + numpy.cos(kvec[2] - kvec[0]))
                    for iorb in range(norb):
                        for jorb in range(norb):
                            if iorb == jorb:
                                print("{0}".format(ek), file=f)  # Real part
                            else:
                                print("0.0", file=f)  # Real part
                    for iorb in range(norb*norb):
                        print("0.0", file=f)  # Imaginary part


def _generate_bethe_lattice_model(nelec, norb, t, nk, f):
    """
    Compute a list of H(k) for bethe lattice model with t

    Parameters
    ----------
    nelec : int
            Number of electrons
    norb : int
            Number of orbitals
    t : float
            Nearest neighbor hopping
    nk : int
            Number of k divisions along each axis
    f : file
            Output file
    """

    nkbz = nk

    print("\n    Total number of k =", str(nkbz))

    #
    # Write General-Hk formatted file
    #
    print(nkbz, file=f)
    print(nelec, file=f)
    print("1", file=f)
    print("0 0 {0} {1}".format(0, norb), file=f)
    print("1", file=f)
    print("0 0 {0} {1} 0 0".format(0, norb), file=f)
    print("1 {0}".format(norb), file=f)

    #
    # If Bethe lattice, set k-weight manually to generate semi-circular DOS
    #
    for i0 in range(nk):
        ek = float(2*i0 + 1 - nk) / float(nk)
        wk = numpy.sqrt(1.0 - ek**2)
        print("{0}".format(wk), file=f)
    for i0 in range(nk):
        ek = 2.0 * t * float(2*i0 + 1 - nk) / float(nk)
        for iorb in range(norb):
            for jorb in range(norb):
                if iorb == jorb:
                    print("{0}".format(ek), file=f)  # Real part
                else:
                    print("0.0", file=f)  # Real part
        for iorb in range(norb*norb):
            print("0.0", file=f)  # Imaginary part


def _generate_wannier90_model(nkdiv, params, f):
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


def para2noncol(params):
    """
    Duplicate orbitals when we perform non-colinear DMFT from the colinear DFT calculation.
    """
    ncor = params["model"]['ncor']
    norb_list = re.findall(r'\d+', params["model"]["norb"])
    norb = [int(norb_list[icor]) / 2 for icor in range(ncor)]
    #
    # Read Wannier90 file
    #
    w90 = Wannier90Converter(seedname=params["model"]["seedname"])
    nr, rvec, rdeg, nwan, hamr = w90.read_wannier90hr(params["model"]["seedname"] + "_col_hr.dat")
    #
    # Non correlated shell
    #
    ncor += 1
    norb_tot = sum(norb)
    norb.append(nwan - norb_tot)
    #
    # Output new wannier90 file
    #
    with open(params["model"]["seedname"] + "_hr.dat", 'w') as f:

        print("Converted from Para to Non-collinear", file=f)
        print(nwan*2, file=f)
        print(nr, file=f)
        for ir in range(nr):
            print("%5d" % rdeg[ir], end="", file=f)
            if ir % 15 == 14:
                print("", file=f)
        if nr % 15 != 0:
            print("", file=f)
        #
        # Hopping
        #
        for ir in range(nr):
            iorb1 = 0
            iorb2 = 0
            for icor in range(ncor):
                for ispin in range(2):
                    iorb1 -= ispin * norb[icor]
                    for iorb in range(norb[icor]):
                        iorb1 += 1
                        iorb2 += 1
                        jorb1 = 0
                        jorb2 = 0
                        for jcor in range(ncor):
                            for jspin in range(2):
                                jorb1 -= jspin * norb[jcor]
                                for jorb in range(norb[jcor]):
                                    jorb1 += 1
                                    jorb2 += 1
                                    if ispin == jspin:
                                        hamr2 = hamr[ir][jorb1-1, iorb1-1]
                                    else:
                                        hamr2 = 0.0 + 0.0j
                                    print("%5d%5d%5d%5d%5d%12.6f%12.6f" %
                                          (rvec[ir, 0], rvec[ir, 1], rvec[ir, 2], jorb2, iorb2,
                                           hamr2.real, hamr2.imag), file=f)


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

    def generate_model_file(self):
        pass


class BetheModel(LatticeModel):
    """
    Bethe-lattice model
    """
    def __init__(self, params):
        super(BetheModel, self).__init__(params)
        self._nk = params["system"]["nk"]

    @classmethod
    def name(cls):
        return 'bethe'

    @classmethod
    def spatial_dim(cls):
        return 1

    def nkdiv(self):
        return (self._nk, 1, 1)

    def generate_model_file(self):
        from pytriqs.applications.dft.converters.hk_converter import HkConverter

        p = self._params
        seedname = p["model"]["seedname"]
        with open(seedname+'.inp', 'w') as f:
            _generate_bethe_lattice_model(p['model']['nelec'], int(p['model']['norb']), p['model']['t'], p['system']['nk'], f)

        # Convert General-Hk to SumDFT-HDF5 format
        converter = HkConverter(filename=seedname + ".inp", hdf_filename=seedname+".h5")
        converter.convert_dft_input(weights_in_file=True)
        os.remove(seedname + ".inp")


class NNNHoppingModel(LatticeModel):
    """
    Nearest-neighbor and next-nearest neighbor hopping model
    """

    def __init__(self, params):
        super(NNNHoppingModel, self).__init__(params)

        self._nk = params["system"]["nk"]

    @classmethod
    def name(cls):
        if cls.spatial_dim() == 1:
            return 'chain'
        elif cls.spatial_dim() == 2:
            return 'square'
        elif cls.spatial_dim() == 3:
            return 'cubic'

    def nkdiv(self):
        return (self._nk,) * self.__class__.spatial_dim() + (1,) * (3-self.__class__.spatial_dim())

    @classmethod
    def spatial_dim(cls):
        raise RuntimeError("spatial_dim must be inherited in a subclass.")

    def generate_model_file(self):
        from pytriqs.applications.dft.converters.hk_converter import HkConverter

        p = self._params
        seedname = p["model"]["seedname"]
        with open(seedname+'.inp', 'w') as f:
            _generate_nnn_lattice_model(self.__class__.spatial_dim(),
                p['model']['nelec'], int(p['model']['norb']), p['model']['t'], p['model']["t'"], p['system']['nk'], f)

        # Convert General-Hk to SumDFT-HDF5 format
        converter = HkConverter(filename=seedname + ".inp", hdf_filename=seedname+".h5")
        converter.convert_dft_input(weights_in_file=False)
        os.remove(seedname + ".inp")


class ChainModel(NNNHoppingModel):
    def __init__(self, params):
        super(ChainModel, self).__init__(params)

    @classmethod
    def spatial_dim(cls):
        return 1


class SquareModel(NNNHoppingModel):
    def __init__(self, params):
        super(SquareModel, self).__init__(params)

    @classmethod
    def spatial_dim(cls):
        return 2


class CubicModel(NNNHoppingModel):
    def __init__(self, params):
        super(CubicModel, self).__init__(params)

    @classmethod
    def spatial_dim(cls):
        return 3


def set_nk(nk, nk0, nk1, nk2):
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


class Wannier90Model(NNNHoppingModel):
    def __init__(self, params):
        super(Wannier90Model, self).__init__(params)

        self._nkdiv = set_nk(params["system"]["nk"],
                             params["system"]["nk0"],
                             params["system"]["nk1"],
                             params["system"]["nk2"])

    @classmethod
    def name(self):
        return 'wannier90'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        from pytriqs.applications.dft.converters.wannier90_converter import Wannier90Converter

        #
        # non_colinear flag is used only for the case that COLINEAR DFT calculation
        #
        p = self._params
        seedname = self._params["model"]["seedname"]
        if p["model"]["spin_orbit"]:
            p["model"]["non_colinear"] = False
        #
        if p["model"]["non_colinear"]:
            para2noncol(p)
        with open(seedname+'.inp', 'w') as f:
            _generate_wannier90_model(self.nkdiv(), p, f)
        # Convert General-Hk to SumDFT-HDF5 format
        converter = Wannier90Converter(seedname=seedname)
        converter.convert_dft_input()


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
        self._nkdiv = set_nk(params["system"]["nk"],
                             params["system"]["nk0"],
                             params["system"]["nk1"],
                             params["system"]["nk2"])

    @classmethod
    def name(self):
        return 'external'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        pass


all_lattice_models = [ChainModel, SquareModel, CubicModel, BetheModel, Wannier90Model, ExternalModel]


