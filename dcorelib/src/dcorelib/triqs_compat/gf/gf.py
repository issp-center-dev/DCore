from copy import deepcopy
import numpy as np
import h5py
import operator

from dcorelib.sparse_gf.basis import matsubara_sampling, tau_sampling, finite_temp_basis

from .meshes import MeshImFreq, MeshImTime, MeshLegendre, MeshIR, MeshReFreq
from ..h5.archive import register_class
from ..plot.protocol import clip_array
from . import plot
from . import meshes

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
            mesh_type = MeshImFreq, target_shape=None):
            self.name = name

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
            if indices is not None and not isinstance(indices, GfIndices):
                raise ValueError("Invalid indices!"+str(indices))
            # At this point, indices is None or an object of GfIndices

            # Determine mesh_type
            if mesh is not None:
                mesh_type = type(mesh)

            # Determine the number of points for freq/time
            if n_points is None:
                if mesh is not None:
                    n_points = mesh.points.size//mesh_type.n_points_fact()
                elif data is not None:
                    n_points = data.shape[0]//mesh_type.n_points_fact()
                elif mesh_type.default_n_points() > 0:
                    n_points = mesh_type.default_n_points()
                else:
                    raise RuntimeError("Failed to determine n_points!")
            
            # Then, dertermine target_shape
            if target_shape is None:
                if indices is not None:
                    target_shape = (len(indices[0]), len(indices[1]))
                elif data is not None:
                    target_shape = data.shape[1:]
                else:
                    raise RuntimeError("Failed to determine target_shape!")
            self.target_shape = tuple(target_shape)

            # beta
            if beta is None:
                if mesh is not None:
                    beta = mesh.beta
                else:
                    raise RuntimeError("Failed to determine beta!")
            self._beta = beta
            
            if statistic is None:
                if mesh is not None:
                    statistic = mesh.statistic
                else:
                    raise RuntimeError("Failed to determine statistic!")
            self._statistic = statistic
            
            # At this point, all necessary information must be set.

            # Construct mesh
            if mesh is None:
                if isinstance(mesh_type, GfImFreq):
                    mesh = mesh_type(beta, statistic=statistic, n_points=n_points//2)
                else:
                    mesh = mesh_type(beta, statistic=statistic, n_points=n_points)
            self.mesh = mesh
            
            # Construct data
            if data is None:
                # Try to figure the shape of data for indices
                data = np.zeros((n_points * mesh_type.n_points_fact(),) + target_shape, dtype=np.complex128)
            self.data = data


            # Construct indices
            if indices is None:
                left_indices = list(map(str, np.arange(self.target_shape[0])))
                right_indices = list(map(str, np.arange(self.target_shape[1])))
                indices = GfIndices([left_indices, right_indices])
            self.indices = indices

        delegate(self, **kw)
    
    def zero(self):
        """ Fill data with zero """
        self.data[...] = 0
    
    def copy(self):
        """ Return a deep copy of self """
        return deepcopy(self)
    
    def to_triqs(self):
        """ Transform to a TRIQS Gf object (data is copied) """
        return NotImplemented

    @classmethod
    def default_n_points():
        return -1

    @classmethod
    def from_triqs(cls, triqs_gf):
        """ Transform a TRIQS Gf object (data is copied) to an instance of the corresponding Gf subclass """
        from triqs.gf import meshes as tmeshes
        cls_ = {
            tmeshes.MeshImFreq: GfImFreq,
            tmeshes.MeshReFreq: GfReFreq,
            tmeshes.MeshLegendre: GfLegendre
        }[type(triqs_gf.mesh)]
        return cls_.from_triqs(triqs_gf)


    def __lshift__(self, A):
        """ Substitute a new gf object (copy) """
        if isinstance(A, Gf):
            for name in ['target_shape', 'name', '_beta', '_statistic']:
                self.__setattr__(name, A.__getattribute__(name))
            self.data[...] = A.data
        elif isinstance(A, np.ndarray):
            if A.ndim == 3:
                self.data[...] = A
            elif A.ndim == 2:
                self.data[...] = A[None,:,:]
            else:
                raise RuntimeError("Invalid ndarray A!")
        elif isinstance(A, LazyExpression):
            self.data[...] = A.evaluate(self).data
        else:
            raise RuntimeError(f"Invalid type of A! {type(A)}")
    
    @property
    def shape(self):
        return self.data.shape

    def __getitem__(self, idx):
        assert isinstance(idx, tuple) and len(idx) == 2
        data_view = self.data[:, idx[0], idx[1]].reshape((-1,1,1))
        g_view = Gf(beta=self._beta, statistic=self._statistic,
            mesh=self.mesh, data=data_view)
        return g_view
    
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

    def invert(self):
        self.data[...] = np.linalg.inv(self.data)

    def conjugate(self):
        g = self.copy()
        g.data[...] = g.data.conjugate()
        return g
    
    def transpose(self):
        g = self.copy()
        g.data[...] = g.data.transpose((0,2,1))
        return g

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
        return gf_subclasses[type(dict['mesh'])](
            data = dict['data'],
            mesh = dict['mesh'],
            beta = dict['mesh'].beta,
            #statistic = dict['mesh'].statistic,
        )

    def __iadd__(self, other):
        self.data[...] += other.data
        return self    

    def __add__(self, other):
        return self.__add_sub__(other, operator.iadd)

    def __sub__(self, other):
        return self.__add_sub__(other, operator.isub)

    def __add_sub__(self, other, op):
        res = self.copy()
        if type(self) == type(other):
            op(res.data, other.data)
        elif isinstance(other, np.ndarray):
            if other.ndim == 3:
                op(res.data, other.data)
            elif other.ndim == 2:
                op(res.data, other.data[None,:,:])
            else:
                raise RuntimeError("Invalid ndarray!")
        return res

    def __isub__(self, other):
        data_ = other.data if isinstance(other, Gf) else other
        self.data -= data_
        return self
    
    def __mul__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        res = self.copy()
        res.data *= other
        return res
    
    __rmul__ = __mul__

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

    #-----------------------------plot protocol -----------------------------------

    def _plot_(self, opt_dict):
        """ Implement the plot protocol"""
        return plot.dispatcher(self)(self, opt_dict)

    def x_data_view(self, x_window=None, flatten_y=False):
        """Helper method for getting a view of the data.

        Parameters
        ----------

        x_window : optional
            The window of x variable (omega/omega_n/t/tau) for which data is requested.
        flatten_y: bool, optional
            If the Greens function is of size (1, 1) flatten the array as a 1d array.

        Returns
        -------

        (X, data) : tuple
            X is a 1d numpy array of the x variable inside the window requested.
            data is a 3d numpy array of dim (:,:, len(X)), the corresponding slice of data.
            If flatten_y is True and dim is (1, 1, *) it returns a 1d numpy array.
        """

        X = [x.imag for x in self.mesh] if isinstance(self.mesh, meshes.MeshImFreq) \
            else [x for x in self.mesh]

        X, data = np.array(X), self.data
        if x_window:
            # the slice due to clip option x_window
            sl = clip_array(X, *x_window) if x_window else slice(len(X))
            X, data = X[sl],  data[sl, :, :]
        if flatten_y and data.shape[1:3] == (1, 1):
            data = data[:, 0, 0]
        return X, data


class GfImFreq(Gf):
    #def __init__(self, **kw):
        #if 'n_points' not in kw:
            #kw['n_points'] = 1025
        #super().__init__(**kw)

    def __lshift__(self, g):
        """Set from GfIR instance"""
        if not isinstance(g, GfIR):
            return super().__lshift__(g)
        smpl = matsubara_sampling(g.basis, sampling_points=self.mesh.points)
        self.data[...] = smpl.evaluate(g.data, axis=0)

    def to_triqs(self):
        from triqs.gf import GfImFreq as _GfImFreq
        return _GfImFreq(
            beta=self.mesh.beta,
            statistic=self.mesh.statistic,
            n_points=self.data.shape[0]//2,
            data=self.data.copy()
            )
    
    @classmethod
    def from_triqs(cls, other):
        return cls(
            beta=other.mesh.beta,
            statistic=other.mesh.statistic,
            n_points=other.data.shape[0]//2,
            data=other.data.copy()
        )

    def inverse(self):
        inv_g = self.copy()
        inv_g.data[...] = np.linalg.inv(self.data)
        return inv_g
    
    def density(self, basis=None):
        if basis is None:
            basis = finite_temp_basis(self._beta, self._statistic)
        smpl = matsubara_sampling(basis, sampling_points=self.mesh.points)
        gl = smpl.fit(self.data, axis=0)
        gbeta = tau_sampling(basis, sampling_points=[self._beta]).evaluate(gl, axis=0)
        return -gbeta[0,:,:].T

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
    

class GfReFreq(Gf):
    def __init__(self, **kw): # enforce keyword only policy 
        if 'window' in kw.keys():
            kw['mesh'] = MeshReFreq(kw['window'][0], kw['window'][1], kw['n_points'])
            del kw['window']
        super().__init__(**kw)

    def to_triqs(self):
        from triqs.gf import GfReFreq as _Gf
        from triqs.gf.meshes import MeshReFreq as _Mesh
        return _Gf(
            data=self.data.copy(),
            mesh = _Mesh(
                self.mesh.points[0], self.mesh.points[-1],
                self.mesh.points.size)
            )
    
    @classmethod
    def from_triqs(cls, other):
        points = np.array([p for p in other.mesh])
        return cls(
            data=other.data.copy(),
            mesh = MeshReFreq(points[0], points[-1], points.size))

    def set_from_pade(self, gm, n_points, freq_offset=0.0):
        """ Set values using Pade approximant """
        assert isinstance(gm, GfImFreq)
        from ..utility.pade_approximants import PadeApproximant
        N = gm.data.shape[1]
        nw = gm.data.shape[0]//2
        idx = range(nw-n_points, nw+n_points)
        z = gm.mesh.values()[idx]
        for i in range(N):
            for j in range(N):
                pade = PadeApproximant(z, gm.data[idx,i,j])
                self.data[:,i,j] = pade(self.mesh.values()+1j*freq_offset)


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

gf_subclasses = {
    MeshImFreq: GfImFreq,
    MeshReFreq: GfReFreq,
    MeshImTime: GfImTime,
    MeshLegendre: GfLegendre,
    MeshIR: GfIR,
}

class LazyExpression(object):
    def __init__(self):
       pass

    def evaluate(g):
        assert np.isinstance(g, GfImFreq)
        return NotImplemented


class LinearExpression(LazyExpression):
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
        elif isinstance(other, GfImFreq):
            val = self.evaluate(other)
            res = other.copy()
            res.data[...] = np.einsum('wij,wjk->wik', val.data, other.data, optimize=True)
            return res
        else:
            return NotImplemented

    def __add__(self, other):
        return self.__add_sub__(other, operator.add)
    
    __radd__ = __add__

    def __sub__(self, other):
        return self.__add_sub__(other, operator.sub)

    def __add_sub__(self, other, op):
        if np.isscalar(other):
            assert np.isscalar(self._a0)
            return LinearExpression(op(self._a0, other), self._a1)
        elif isinstance(other, np.ndarray):
            a0_ = self._a0 if isinstance(self._a0, np.ndarray) else self._a0 * np.identity(other.shape[0])
            return LinearExpression(op(a0_, other), self._a1)
        elif isinstance(other, GfImFreq):
            return op(self.evaluate(other), other)
        else:
            return NotImplemented


    def evaluate(self, g):
        res = g.copy()
        res.zero()

        a0_ = _convert_to_matrix(self._a0, g)
        a1_ = _convert_to_matrix(self._a1, g)

        res.data[...] += a0_[None,:,:]
        w = g.mesh.values()
        res.data[...] += w[:,None,None] * a1_[None,:,:]
        return res

    def inverse(self):
        return InverseLinearExpression(self)

def _convert_to_matrix(a, g):
    if np.isscalar(a):
        nf = g.data.shape[1]
        return a * np.identity(nf)
    else:
        return a

class InverseLinearExpression(LazyExpression):
    """ Inverse of Linear Expression in frequency
    """
    def __init__(self, lin_exp):
        super().__init__()
        assert isinstance(lin_exp, LinearExpression)
        self._lin_exp = lin_exp

    def evaluate(self, g):
        return self._lin_exp.evaluate(g).inverse()
    
    
    def inverse(self):
        return self._lin_exp

# Evalaute to iv (0 + 1*z)
iOmega_n = LinearExpression(0., 1.)
Omega = LinearExpression(0., 1.)

Identity = LinearExpression(1., 0.)


class SpectralModel(LazyExpression):
    """
    Diagonal Green's function generated by a model spectrum
    """
    def __init__(self, omega_min, omega_max, rho_omega):
        """
        Args:
            omega_min (float): lower bound for the spectral function
            omega_max (float): upper bound for the spectral function
            rho_omega (function): Return the value of the spectral function
        """
        super().__init__()
        self._omega_min = omega_min
        self._omega_max = omega_max
        self._rho_omega = rho_omega

    def evaluate(self, g):
        wmax = max(abs(self._omega_max), abs(self._omega_min))
        beta = g.mesh.beta
        nf = g.data.shape[1]
        lambda_ = beta * wmax
        basis = finite_temp_basis(beta, g.mesh.statistic, lambda_, eps=1e-7)
        gl = -basis.v.overlap(self._rho_omega) * basis.s
        giv = matsubara_sampling(basis, sampling_points=g.mesh.points).evaluate(gl)
        return np.einsum('w,ij->wij', giv, np.identity(nf), optimize=True)


class SemiCircular(SpectralModel):
    """
    Diagonal Green's function generated by a model spectrum
    A(omega) = 2\sqrt(D^2-omega^2)/(pi D^2)
    """
    def __init__(self, D, coeff=1):
        """
        Args:
            D (float): Half band width
        """
        super().__init__(-D, D, lambda x: coeff*2*np.sqrt(D**2-x**2)/(np.pi*D**2))
        self.D = D
        self.coeff = coeff
    
    def __mul__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        return SemiCircular(self.D, other*self.coeff)

    __rmul__ = __mul__