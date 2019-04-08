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
import sys
import numpy
import scipy
from itertools import product
from pytriqs.archive.hdf_archive import HDFArchive

from .tools import pauli_matrix

def create_lattice_model(params):
    model_name =  params["model"]["lattice"]

    for c in all_lattice_models:
        if model_name == c.name():
            return c(params)

    raise RuntimeError('Unknown lattice model name: ' + model_name)


def _generate_bethe_lattice_model(norb, t, nk):
    """
    Compute a list of H(k) and weights for bethe lattice model with t

    Parameters
    ----------
    norb : int
            Number of orbitals
    t : float
            Nearest neighbor hopping
    nk : int
            Number of k divisions along each axis
    """

    #
    # If Bethe lattice, set k-weight manually to generate semi-circular DOS
    #
    weight = numpy.empty((nk,))
    for ik in range(nk):
        ek = float(2*ik + 1 - nk) / float(nk)
        weight[ik] = numpy.sqrt(1.0 - ek**2)

    Hk = numpy.zeros((nk, norb, norb), dtype=complex)
    for ik in range(nk):
        ek = 2.0 * t * float(2*ik + 1 - nk) / float(nk)
        for iorb in range(norb):
            Hk[ik, iorb, iorb] = ek

    return Hk, weight


def _call_Hk_converter(seedname, nelec, norb, Hk, weight):
    """
    Compute a list of H(k) for bethe lattice model with t

    Parameters
    ----------
    nelec : int
            Number of electrons
    norb : int
            Number of orbitals
    Hk : [nkbz, :, :]
            H(k)
    weight : [nkbz] or None
            weight for k
    """
    from converters.hk_converter import HkConverter

    nkbz = Hk.shape[0]
    assert nelec <= norb
    assert Hk.shape == (nkbz, norb, norb)
    assert weight is None or weight.shape == (nkbz,)

    use_weight = not weight is None

    print("\n    Total number of k =", str(nkbz))

    with open(seedname+'.inp', 'w') as f:
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

        # Write weight
        if use_weight:
            for ik in range(nkbz):
                print("{0}".format(weight[ik]), file=f)

        # Write H(k)
        for ik in range(nkbz):
            for iorb, jorb in product(range(norb), repeat=2):
                print("{}".format(Hk[ik,iorb,jorb].real), file=f)
            for iorb, jorb in product(range(norb), repeat=2):
                print("{}".format(Hk[ik,iorb,jorb].imag), file=f)

    converter = HkConverter(filename=seedname + ".inp", hdf_filename=seedname+".h5")
    converter.convert_dft_input(weights_in_file=use_weight)
    os.remove(seedname + ".inp")


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

    def Hk(self, kvec):
        raise RuntimeError("Hk is ill-defied for BetheModel")

    def generate_model_file(self):
        p = self._params
        seedname = p["model"]["seedname"]
        Hk, weight = _generate_bethe_lattice_model(int(p['model']['norb']), p['model']['t'], p['system']['nk'])
        _call_Hk_converter(seedname, p['model']['nelec'], int(p['model']['norb']), Hk, weight)

        if p['model']['spin_orbit']:
            from .manip_database import turn_on_spin_orbit
            print('')
            print('Turning on spin_orbit...')
            turn_on_spin_orbit(seedname + '.h5', seedname + '.h5')


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

    def Hk(self, kvec):
        spatial_dim = self.__class__.spatial_dim()
        t = self._params['model']['t']
        tp = self._params['model']["t'"]
        norb = int(self._params['model']['norb'])

        Hk = numpy.zeros((norb, norb), dtype=complex)
        if spatial_dim == 1:
            ek = 2.0*t*numpy.cos(kvec[0]) + 2*tp*numpy.cos(2.0*kvec[0])
            for iorb in range(norb):
                Hk[iorb, iorb] = ek
        elif spatial_dim == 2:
           ek = 2.0*t*(numpy.cos(kvec[0]) + numpy.cos(kvec[1])) \
                + 2.0*tp*(numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1]))
           for iorb in range(norb):
               Hk[iorb, iorb] = ek
        elif spatial_dim == 3:
           ek = 2*t*(numpy.cos(kvec[0]) + numpy.cos(kvec[1]) + numpy.cos(kvec[2])) \
                + 2*tp*(numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1])
                        + numpy.cos(kvec[1] + kvec[2]) + numpy.cos(kvec[1] - kvec[2])
                        + numpy.cos(kvec[2] + kvec[0]) + numpy.cos(kvec[2] - kvec[0]))
           for iorb in range(norb):
               Hk[iorb, iorb] = ek

        if self._params['model']['spin_orbit']:
            return scipy.linalg.block_diag((Hk, Hk))
        else:
            return (Hk, Hk)

    def generate_model_file(self):
        p = self._params
        seedname = p['model']['seedname']
        spin_orbit = p['model']['spin_orbit']
        spatial_dim = self.__class__.spatial_dim()
        norb = int(p['model']['norb'])
        nk = p['system']['nk']
        nkbz = nk**spatial_dim
        Hk = numpy.zeros((nkbz, norb, norb), dtype=complex)
        nkdiv = self.nkdiv()

        # Since Hk_converter does support SO=1, we create a model file for a spinless model.
        if spin_orbit:
            for ik, (ik0, ik1, ik2) in enumerate(product(range(nkdiv[0]), range(nkdiv[1]), range(nkdiv[2]))):
                Hk_tmp = self.Hk( 2*numpy.pi*numpy.array((ik0, ik1, ik2))/float(nk))
                assert numpy.allclose(Hk_tmp[0:norb, 0:norb], Hk_tmp[norb:2*norb, norb:2*norb])
                Hk[ik, :, :] = Hk_tmp[0:norb, 0:norb]
        else:
            for ik, (ik0, ik1, ik2) in enumerate(product(range(nkdiv[0]), range(nkdiv[1]), range(nkdiv[2]))):
                Hk_tmp = self.Hk( 2*numpy.pi*numpy.array((ik0, ik1, ik2))/float(nk))
                assert numpy.allclose(Hk_tmp[0], Hk_tmp[1])
                Hk[ik, :, :] = Hk_tmp[0]
        _call_Hk_converter(seedname, p['model']['nelec'], int(p['model']['norb']), Hk, None)

        if p['model']['spin_orbit']:
            from .manip_database import turn_on_spin_orbit
            print('')
            print('Turning on spin_orbit...')
            turn_on_spin_orbit(seedname + '.h5', seedname + '.h5')


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


class Wannier90Model(NNNHoppingModel):
    def __init__(self, params):
        super(Wannier90Model, self).__init__(params)

        self._nkdiv = _set_nk(params["system"]["nk"],
                             params["system"]["nk0"],
                             params["system"]["nk1"],
                             params["system"]["nk2"])

    @classmethod
    def name(self):
        return 'wannier90'

    def nkdiv(self):
        return self._nkdiv

    def generate_model_file(self):
        from converters.wannier90_converter import Wannier90Converter

        #
        # non_colinear flag is used only for the case that COLINEAR DFT calculation
        #
        p = self._params
        seedname = self._params["model"]["seedname"]
        if p["model"]["spin_orbit"] and p["model"]["non_colinear"]:
            raise RuntimeError("The options spin_orbit and non_colinear are not compatible!")
        with open(seedname+'.inp', 'w') as f:
            _generate_wannier90_model(self.nkdiv(), p, f)

        # Convert General-Hk to SumDFT-HDF5 format
        converter = Wannier90Converter(seedname=seedname)
        converter.convert_dft_input()

        if p["model"]["spin_orbit"] or p["model"]["non_colinear"]:
            if p['model']['non_colinear']:
                from .manip_database import turn_on_spin_orbit
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
        self._nkdiv = _set_nk(params["system"]["nk"],
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

def _drop_small_vals(z, eps=1e-10):
    z_real = z.real if numpy.abs(z.real) > eps else 0.0
    z_imag = z.imag if numpy.abs(z.imag) > eps else 0.0
    return complex(z_real, z_imag)

def print_spin_orbital_matrix(mat, file, print_offset=0):
    norb = mat.shape[1]

    for isp in range(2):
        for iorb in range(norb):
            print(' '*print_offset, end='')
            for jsp in range(2):
                for jorb in range(norb):
                    z = _drop_small_vals(mat[isp, iorb, jsp, jorb])
                    print('({0:>9.2e},{1:>9.2e})'.format(z.real, z.imag), end=' ', file=file)
                if jsp==0:
                    print('  ', file=file, end='')
            print('', file=file)
        print('', file=file) 

def print_local_fields(h5_file, corr_shell_dims=None, subgrp='dft_input'):
    """
    Print local fields of H(R=0)
    :param h5_file: input file for DFTTools
    :param subgrp:
    """

    with HDFArchive(h5_file, 'r') as f:
        Hk = f[subgrp]['hopping'][()]
        SO = f[subgrp]['SO']
        SP = f[subgrp]['SP']
        bz_weights = f[subgrp]['bz_weights'][()]
        n_corr_sh = f[subgrp]['n_corr_shells']
        dims_corr_sh = numpy.array([f[subgrp]['corr_shells'][ish]['dim'] for ish in range(n_corr_sh)])

    if (SO==1 and SP==0) or (SO==0 and SP==1):
        raise RuntimeError("SO={} and SP={} are not supported by DCore!".format(SO, SP))

    nk = Hk.shape[0]
    spin_block_dim = Hk.shape[2]
    if SO==0:
        Hk_ud = numpy.zeros((nk, 2, spin_block_dim, 2, spin_block_dim), dtype=complex)
        for isp in range(2):
            Hk_ud[:, isp, :, isp, :] = Hk[:, 0, :, :]
        Hk_ud = Hk_ud.reshape((nk, 2*spin_block_dim, 2*spin_block_dim))
        norb = spin_block_dim
        num_spin_orb_corr_sh = 2*dims_corr_sh
    else:
        Hk_ud = Hk.reshape((nk, spin_block_dim, spin_block_dim))
        norb = spin_block_dim//2
        num_spin_orb_corr_sh = dims_corr_sh

    # Hk_ud (nk, 2*num_orb, 2*num_orb)
    # Check only H(R=0)
    H0_ud = numpy.einsum('kij,k->ij', Hk_ud, bz_weights)
    H0_ud = H0_ud.reshape((2, norb, 2, norb))

    print('')
    print('---local fields (w/o local potential)')
    pauli_mat = pauli_matrix()
    for iorb in range(norb):
        # alpha = x, y, z
        h = numpy.empty(3)
        for alpha in range(3):
            h[alpha] = 0.5 * numpy.trace(numpy.dot(H0_ud[:, iorb, :, iorb], pauli_mat[alpha])).real
        print("    orbital {} : hx, hy, hz = {} {} {}".format(iorb, *h))

    print('')
    print('---intra shell Hamiltonian (w/o local potential)')
    # (spin, orb, spin, orb) => (orb, spin, orb, spin)
    H0_ud2 = H0_ud.transpose((1, 0, 3, 2)).reshape((2*norb, 2*norb))
    offset = 0
    print('    ordering: up ... up, down ... down')
    for ish in range(n_corr_sh):
        print(' '*4, 'corr_shell=', ish)
        block_size = num_spin_orb_corr_sh[ish]
        H0_block = H0_ud2[offset:offset+block_size, offset:offset+block_size]
        # (orb, spin, orb, spin) => (spin, orb, spin, orb)
        H0_block = H0_block.reshape((block_size//2, 2, block_size//2, 2)).transpose((1, 0, 3, 2))
        print_spin_orbital_matrix(H0_block, sys.stdout, 6)
        evals, evecs = numpy.linalg.eigh(H0_block.reshape((block_size, block_size)))
        print('      eigenvalues: ', evals)
        print('')
        offset += block_size

all_lattice_models = [ChainModel, SquareModel, CubicModel, BetheModel, Wannier90Model, ExternalModel]
