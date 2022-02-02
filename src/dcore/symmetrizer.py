import numpy
from numpy.core.numeric import identity
from scipy.linalg import expm
from ._dispatcher import BlockGf

from .tools import pauli_matrix, symmetrize, symmetrize_spin

"""
Symmetrizer for local Green's function/self-energy (BlockGf)
"""
class LocalSymmetrizer(object):
    def __init__(self):
        super().__init__()

    def __call__(self, gf):
        """ Symmetrize local Green's function """
        return NotImplemented

class LocalSymmetrizerFromGenerators(LocalSymmetrizer):
    def __init__(self, symm_generators):
        """
        symmetrizers: list of generators of symmetrization
            Each generator is a dict (key: block name, value: a 2D array)
        """
        assert isinstance(symm_generators, list)
        self._symm_gen = symm_generators

    def __call__(self, gf):
        """ Symmetrize local Green's function """
        assert isinstance(gf, BlockGf)
        return symmetrize(gf, self._symm_gen)

class PMSymmSpinOrbitOff(LocalSymmetrizer):
    def __init__(self):
        super().__init__()
    
    def __call__(self, gf):
        """ Symmetrize local Green's function """
        assert isinstance(gf, BlockGf)
        gf_copy = gf.copy()
        symmetrize_spin(gf_copy)
        return gf_copy


def _pm_symm_gen(norb, transverse: bool):
    theta = numpy.pi
    sx, sy, sz = pauli_matrix()
    _gen = lambda s: numpy.einsum('pq,ST->SpTq', numpy.identity(norb), expm(-0.5j*theta*s)).reshape(2*norb, 2*norb)
    if transverse:
        gen_z = _gen(sx)
        gen_x = _gen(sy)
        gen_y = _gen(sz)
        return [{"ud": gen_x}, {"ud": gen_y}, {"ud": gen_z}]
    else:
        gen_z = _gen(sx)
        return [{"ud": gen_z}]


def pm_symmetrizer(norb: int, spin_orbit: bool, transverse: bool):
    if spin_orbit:
        return LocalSymmetrizerFromGenerators(_pm_symm_gen(norb, transverse))
    else:
        return PMSymmSpinOrbitOff()