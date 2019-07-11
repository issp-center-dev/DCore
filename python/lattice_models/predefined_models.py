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
import numpy
import scipy
from itertools import product
from pytriqs.archive.hdf_archive import HDFArchive

from .base import LatticeModel

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
            Number of orbitals (not physical orbitals actually the size of H(k))
    Hk : [nkbz, :, :]
            H(k)
    weight : [nkbz] or None
            weight for k
    """
    from ..converters.hk_converter import HkConverter

    nkbz = Hk.shape[0]
    #assert nelec <= norb
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


class BetheModel(LatticeModel):
    """
    Bethe-lattice model
    """
    def __init__(self, params):
        super(BetheModel, self).__init__(params)
        self._nk = params["model"]["nk"]

    @classmethod
    def name(cls):
        return 'bethe'

    @classmethod
    def spatial_dim(cls):
        return 1

    def nkdiv(self):
        return (self._nk, 1, 1)

    @classmethod
    def is_Hk_supported(cls):
        return False

    def Hk(self, kvec):
        raise RuntimeError("Hk is ill-defied for BetheModel")

    def generate_model_file(self):
        p = self._params
        seedname = p["model"]["seedname"]
        Hk, weight = _generate_bethe_lattice_model(int(p['model']['norb']), p['model']['t'], p['model']['nk'])
        _call_Hk_converter(seedname, p['model']['nelec'], int(p['model']['norb']), Hk, weight)

        if p['model']['spin_orbit']:
            from ..manip_database import turn_on_spin_orbit
            print('')
            print('Turning on spin_orbit...')
            turn_on_spin_orbit(seedname + '.h5', seedname + '.h5')


class NNNHoppingModel(LatticeModel):
    """
    Nearest-neighbor and next-nearest neighbor hopping model
    """

    def __init__(self, params):
        super(NNNHoppingModel, self).__init__(params)

        self._nk = params["model"]["nk"]

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
            return scipy.linalg.block_diag(Hk, Hk)
        else:
            return (Hk, Hk)

    def generate_model_file(self):
        p = self._params
        seedname = p['model']['seedname']
        spin_orbit = p['model']['spin_orbit']
        spatial_dim = self.__class__.spatial_dim()
        norb = int(p['model']['norb'])
        nk = p['model']['nk']
        nkbz = nk**spatial_dim
        Hk = numpy.zeros((nkbz, norb, norb), dtype=complex)
        nkdiv = self.nkdiv()

        # Since Hk_converter does support SO=1, we create a model file for a spinless model.
        if spin_orbit:
            for ik, (ik0, ik1, ik2) in enumerate(product(range(nkdiv[0]), range(nkdiv[1]), range(nkdiv[2]))):
                Hk_ud = self.Hk( 2*numpy.pi*numpy.array((ik0, ik1, ik2))/float(nk))
                assert numpy.allclose(Hk_ud[0:norb, 0:norb], Hk_ud[norb:2*norb, norb:2*norb])
                Hk[ik, :, :] = Hk_ud[0:norb, 0:norb]
        else:
            for ik, (ik0, ik1, ik2) in enumerate(product(range(nkdiv[0]), range(nkdiv[1]), range(nkdiv[2]))):
                Hk_up_down = self.Hk( 2*numpy.pi*numpy.array((ik0, ik1, ik2))/float(nk))
                assert numpy.allclose(Hk_up_down[0], Hk_up_down[1])
                Hk[ik, :, :] = Hk_up_down[0]
        _call_Hk_converter(seedname, p['model']['nelec'], int(p['model']['norb']), Hk, None)

        if p['model']['spin_orbit']:
            from ..manip_database import turn_on_spin_orbit
            print('')
            print('Turning on spin_orbit...')
            turn_on_spin_orbit(seedname + '.h5', seedname + '.h5')


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
        """

        n_k = kvec.shape[0]
        assert kvec.shape[1] == 3

        norb_sh = numpy.array(params['model']['norb_corr_sh'])

        assert len(norb_sh) == 1
        assert params['model']['ncor'] == 1

        #
        # Energy band
        #
        norb = norb_sh[0]
        dim_Hk = 2*norb if params['model']['spin_orbit'] else norb
        n_spin = 1
        n_orbitals = numpy.ones((n_k, n_spin), dtype=int) * dim_Hk
        hopping = numpy.zeros((n_k, n_spin, dim_Hk, dim_Hk), complex)
        if params['model']['spin_orbit']:
            for ik in range(n_k):
                hopping[ik, 0, :, :] = self.Hk(kvec[ik])
        else:
            for ik in range(n_k):
                # Copy only the up component
                hopping[ik, 0, :, :] = self.Hk(kvec[ik])[0]

        #
        # proj_mat is (norb*norb) identities at each correlation shell
        #
        proj_mat = numpy.zeros([n_k, n_spin, 1, dim_Hk, dim_Hk], complex)
        proj_mat[:, :, 0, 0:dim_Hk, 0:dim_Hk] = numpy.identity(dim_Hk, complex)

        #
        # Output them into seedname.h5
        #
        with HDFArchive(params['model']['seedname'] + '.h5', 'a') as f:
            if not (bands_data in f):
                f.create_group(bands_data)
            f[bands_data]['hopping'] = hopping
            f[bands_data]['n_k'] = n_k
            f[bands_data]['n_orbitals'] = n_orbitals
            f[bands_data]['proj_mat'] = proj_mat
        print('    Done')


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

