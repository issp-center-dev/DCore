from copy import deepcopy, copy
import numpy as np

from dcore.backend.sparse_gf.basis import matsubara_sampling, tau_sampling

from .meshes import MeshImFreq, MeshImTime, MeshLegendre, MeshIR
from ..h5.archive import register_class
import h5py

def _to_fixed_length_utf8_array(str_list):
    length = int(np.amax([len(x) for x in str_list]))
    dt = h5py.string_dtype(encoding='utf-8', length=length)
    return np.array(str_list, dtype=dt)

def _to_utf8_strings(str_list):
    return [(s.decode(encoding='utf-8') if isinstance(s, bytes) else s) for s in str_list]
class GfIndices:
    def __init__(self, indices):
        """GfIndices

        Args:
            indices (list): list of list of str
        """
        assert isinstance(indices, list)
        self._data = [_to_utf8_strings(x) for x in indices]
    
    @property
    def data(self):
        return self._data
    
    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        assert len(self._data) == 2
        group.create_group(key)
        group[key]['left'] = _to_fixed_length_utf8_array(self._data[0])
        group[key]['right'] = _to_fixed_length_utf8_array(self._data[1])
        group[key].write_attr('Format', 'GfIndices')
    
    def __getitem__(self, key):
        return self._data[key]
    
    @classmethod
    def __factory_from_dict__(cls, key, dict) :
        return cls([dict['left'], dict['right']])
    
    def __len__(self):
        return len(self._data[0])


def _is_list_of(objs, expected_elem_type):
    """ Return if objs is a list of objects of expected_elem_type """
    if not isinstance(objs, list):
        return False
    for x in objs:
        if not isinstance(x, expected_elem_type):
            return False
    return True

class Gf(object): 
    """
    Parameters (KEYWORD argument ONLY)
    ----------
    beta: float
        Inverse temperature

    statistic: str
        'Fermion' (default) or 'Boson'

    mesh: Mesh object
        If not given, a MeshImFreq object is constructed.

    data: numpy.array, optional
        The data of the Gf.
        The size of the first dimension is the number of points (frequencies, times).

    name: str
        The name of the Green function. For plotting.
    
    indices: WARNING: The Use of string indices is deprecated!
             GfIndices or list of str(int) or list of list of str(int), optional
             Optional string indices for the target space, to allow e.g. ``['eg', 'eg']``.
             list of list of str/int: the list of indices for each dimension.
             list of str/int: all indices are assumed to be the same for all dimensions.
    
    n_points: int
        DEPRECATED:
        Number of points (frequencies/taus/legendere polys).
        For imaginary-frequencies Green's funciton, the first dimension of `data' will be 2*n_points
        because the data includes both of positive and negative frequencies.
        If this option is given, data and mesh must be None.
    """
    def __init__(self, **kw): # enforce keyword only policy 
        def delegate(self, data=None, name='', beta=None, statistic='Fermion', mesh=None, indices=None, n_points=None,
            mesh_type = MeshImFreq):
            # Check indices
            if isinstance(indices, np.ndarray) or isinstance(indices, list):
                indices = list(map(str, indices))
            if indices is None:
                pass
            elif _is_list_of(indices, str):
                # List of strings
                indices = GfIndices(2*[indices])
            elif isinstance(indices, list):
                for x in indices:
                    assert _is_list_of(x, str)
            else:
                raise ValueError("Invalid indices!")
            # At this point, indices is None or an object of GfIndices

            # First determine n_points
            if n_points is not None:
                assert data is None
                assert mesh is None
                mesh = mesh_type(beta, statistic=statistic, n_points=n_points)
            
            if data is None:
                # Try to figure the shape of data for indices
                assert indices is not None
                N1, N2, = len(indices[0]), len(indices[1])
                data = np.empty((mesh._points.size, N1, N2), dtype=np.complex128)

            self.data = data
            self.target_shape = self.data[1:]

            self.name = name
            self.beta = beta
            self.statistic = statistic
            if mesh is None:
                if mesh_type == MeshImFreq:
                    mesh = mesh_type(beta, statistic, self.data.shape[0]//2)
                else:
                    mesh = mesh_type(beta, statistic, self.data.shape[0])
            self.mesh = mesh
            if indices is None:
                left_indices = list(map(str, np.arange(self.data.shape[1])))
                right_indices = list(map(str, np.arange(self.data.shape[2])))
                indices = GfIndices([left_indices, right_indices])
            self.indices = indices

        delegate(self, **kw)
    
    def zero(self):
        """ Fill data with zero """
        self.data[...] = 0
    
    def copy(self):
        """ Return a deep copy of self """
        return deepcopy(self)

    def __lshift__(self, A):
        """ Substitute a new gf object (copy) """
        if isinstance(A, Gf):
            for name in ['data', 'target_shape', 'name', 'beta', 'statistic']:
                self.__setattr__(name, copy(A.__getattribute__(name)))
        elif isinstance(A, np.ndarray):
            self.data[...] = A
        elif type(A) in [LinearExpression, InverseLinearExpression]:
            self.data[...] = A.evaluate(self).data
        else:
            raise RuntimeError(f"Invalid type of A! {type(A)}")
    
    @property
    def shape(self):
        return self.data.shape
    
    def from_L_G_R(self, L, G, R):
        """Matrix transform of the target space of a matrix valued Greens function.

        Sets the current Greens function :math:`g_{ab}` to the matrix transform of :math:`G_{cd}`
        using the left and right transform matrices :math:`L_{ac}` and :math:`R_{db}`.

        .. math::
            g_{ab} = \sum_{cd} L_{ac} G_{cd} R_{db}

        Parameters
        ----------

        L : (a, c) ndarray
            Left side transform matrix.
        G : Gf matrix valued target_shape == (c, d)
            Greens function to transform.
        R : (d, b) ndarray
            Right side transform matrix.
        """
        assert L.ndim == 2, "L needs to be two dimensional"
        assert R.ndim == 2, "R needs to be two dimensional"
        assert L.shape[0] == self.data.shape[1], "Dimensions of L and self are not compatible"
        assert L.shape[1] == G.shape[1], "Dimensions of L and G are not compatible"
        assert G.shape[2] == R.shape[0], "Dimensions of G and R are not compatible"
        assert R.shape[1] == self.shape[2], "Dimensions of R and self are not compatible"
        assert G.shape[0] == self.shape[0], "The leading dimensions of G and self are not compatible"

        self.data[...] = np.einsum('ac,wcd,db->wab', L, G.data, R, optimize=True)


    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        group.create_group(key)
        group[key]['data'] = self.data
        group[key]['mesh'] = self.mesh
        group[key]['indices'] = self.indices
        group[key].write_attr('Format', 'Gf')
        assert self.data.ndim == 3

    @classmethod
    def __factory_from_dict__(cls, key, dict) :
        return cls(
            data = dict['data'],
            mesh = dict['mesh'],
            beta = dict['mesh'].beta,
            statistic = dict['mesh'].statistic,
        )
        #return cls([dict['left'], dict['right']])

    def __iadd__(self, other):
        self.data[...] += other.data
        return self    

    def __add__(self, other):
        res = self.copy()
        return res + other

    def __mul__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        res = self.copy()
        res.data *= other
        return res

    def __truediv__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        res = self.copy()
        res /= other
        return res

    def __itruediv__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        self.data /= other
        return self


class GfImFreq(Gf):
    def __lshift__(self, g):
        """Set from GfIR instance"""
        if not isinstance(g, GfIR):
            return super().__lshift__(g)
        smpl = matsubara_sampling(g.basis, sampling_points=self.mesh.points)
        self.data[...] = smpl.evaluate(g.data, axis=0)
    
    def inverse(self):
        inv_g = self.copy()
        inv_g.data[...] = np.linalg.inv(self.data)
        return inv_g

class GfImTime(Gf):
    def __init__(self, **kw): # enforce keyword only policy 
        if 'mesh' not in kw.keys() or kw['mesh'] is None:
            mesh_type = MeshImTime
        super().__init__(**kw, mesh_type=mesh_type)

    def __lshift__(self, g):
        """Set from GfIR instance"""
        if not isinstance(g, GfIR):
            return super().__lshift__(g)
        smpl = tau_sampling(g.basis, sampling_points=self.mesh.points)
        self.data[...] = smpl.evaluate(g.data, axis=0)

class GfLegendre(Gf):
    def __init__(self, data=None, indices=None, beta=None, n_points=None, name=""):
        super().__init__(data=data, indices=indices, beta=beta, mesh=MeshLegendre(n_points), name=name)


class GfIR(Gf):
    def __init__(self, data, basis, beta=None, name=""):
        super().__init__(data=data, indices=None, beta=beta, mesh=MeshIR(basis), name=name)
        self.basis = basis

register_class(GfIndices)
register_class(Gf)
register_class(GfImFreq)
register_class(GfImTime)
register_class(GfLegendre)
register_class(GfIR)

class LazyExpression(object):
    def __init__():
       pass

    def evaluate(g):
        assert np.isinstance(g, GfImFreq)
        return NotImplemented

class LinearExpression(object):
    """Linear Expression in frequency

    a_0 + a_1 * z,
        where z is a frequency.

    a_i is a scalar or a matrix.
    A scalar is interpreted as `scalar * identity matrix`.
    """
    def __init__(self, a0=0., a1=1.):
        super().__init__()
        self._a0 = a0
        self._a1 = a1

    def copy(self):
        return LinearExpression(self._a0, self._a1)

    def __mul__(self, other):
        if np.isscalar(other):
            return LinearExpression(other*self._a0, other*self._a1)
        return NotImplemented

    def __add__(self, other):
        if np.isscalar(other):
            assert np.isscalar(self._a0)
            return LinearExpression(self._a0 + other, self._a1)
        elif isinstance(other, np.ndarray):
            a0_ = self._a0 if isinstance(self._a0, np.ndarray) else self._a0 * np.identity(other.shape[0])
            return LinearExpression(a0_ + other, self._a1)
        elif np.isinstance(other, GfImFreq):
            return self.evaluate() + other
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (-1) * other

    def evaluate(self, g):
        res = g.copy()
        res.zero()

        res.data[...] += self._a0[None,:,:]
        if np.isscalar(self._a1):
            nf = g.data.shape[1]
            a1_ = self._a1 * np.identity(nf)
        else:
            a1_ = self._a1
        iv = 1J*g.mesh.points * np.pi/g.beta
        res.data[...] += iv[:,None,None] * a1_[None,:,:]
        return res

    def inverse(self):
        return InverseLinearExpression(self)

class InverseLinearExpression(object):
    """ Inverse of Linear Expression in frequency
    """
    def __init__(self, lin_exp):
        super().__init__()
        assert isinstance(lin_exp, LinearExpression)
        self._lin_exp = lin_exp

    def evaluate(self, g):
        return inverse(self._lin_exp.evaluate(g))
    
    def inverse(self):
        return self._lin_exp

# Evalaute to iv (0 + 1*z)
iOmega_n = LinearExpression(0., 1.)


def inverse(g):
    """
    Compute inverse of imaginary-frequency Green's function
    """
    return g.inverse()