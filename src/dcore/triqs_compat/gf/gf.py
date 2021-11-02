import numpy as np
import warnings

from . import meshes

class GfIndices:
    def __init__(self, l):
        """

        Args:
            l (list of str): indices
        """
        assert isinstance(l, list)
        for x in l:
            isinstance(x, str)
        self._l = l
    
    def __iter__(self):
        yield from self._l
    
    #def transpose(self):
        #return GfIndices([])

class Gf:
    r""" Greens function container class

    Parameters
    ----------

    mesh: Greens function mesh
          One of the meshes of the module 'meshes'.

    data: numpy.array, optional
          The data of the Greens function.
          Must be of dimension ``mesh.rank + target_rank``.

    target_shape: list of int, optional
                  Shape of the target space.

    is_real: bool
             Is the Greens function real valued?
             If true, and target_shape is set, the data will be real.
             Mutually exclusive with argument ``data``.

    indices: 
             GfIndices or list of str
   
    name: str 
          The name of the Greens function for plotting.

    Notes
    -----

    One of ``target_shape`` or ``data`` must be set, and the other must be `None`. If passing ``data``
    and ``indices`` the ``data.shape`` needs to be compatible with the shape of ``indices``.

    """
    
    _hdf5_data_scheme_ = 'Gf'

    def __init__(self, **kw): # enforce keyword only policy 
        
        #print "Gf construct args", kw

        def delegate(self, mesh, data=None, target_shape=None, indices = None, name = '', is_real = False):
            """
            target_shape and data  : must provide exactly one of them
            """
            self.name = name

            # input check
            assert (target_shape is None) or (data is None), "data and target_shape : one must be None"
            assert (data is None) or (is_real is False), "is_real can not be True if data is not None"
            if target_shape : 
                for i in target_shape : 
                    assert i>0, "Target shape elements must be >0"
     
            # mesh
            assert isinstance(mesh, meshes.all_meshes), "Mesh is unknown. Possible type of meshes are %s" % ', '.join([m.__name__ for m in all_meshes])
            self._mesh = mesh

            # indices
            assert isinstance(indices, (type(None), list)), "Type of indices incorrect : should be None, list of str, or list of list of str"
            if isinstance(indices, list):
                indices = [x if isinstance(x, str) else str(x) for x in indices]
            self._indices = indices # now indices are None or Gfindices 

            # data
            if data is None:
                # if no data, we get the target_shape. If necessary, we find it from of the list of indices
                if target_shape is None : 
                    assert indices, "Without data, target_shape, I need the indices to compute the shape !"
                    target_shape = [ len(x) for x in indices.data]
                # we now allocate the data
                l = mesh.size
                data = np.zeros(list(l) + list(target_shape), dtype = np.float64 if is_real else np.complex128)
            else:
                l = (len(mesh),)
                assert l == data.shape[0:len(l)], "Mismatch between data shape %s and sizes of mesh(es) %s\n " % (data.shape, l)

            # Now we have the data at correct size. Set up a few short cuts
            self._data = data
            len_data_shape = len(self._data.shape) 
            self._target_rank = len_data_shape - 1
            self._rank = len_data_shape - self._target_rank 
            assert self._rank >= 0

            # target_shape. Ensure it is correct in any case.
            assert target_shape is None or tuple(target_shape) == self._data.shape[self._rank:] # Debug only
            self._target_shape = self._data.shape[self._rank:]

            # If no indices was given, build the default ones
            if self._indices is None: 
                self._indices = GfIndices([list(str(i) for i in range(n)) for n in self._target_shape])

            # Check that indices  have the right size
            if self._indices is not None: 
                d,i =  self._data.shape[self._rank:], tuple(len(x) for x in self._indices.data)
                assert (d == i), "Indices are of incorrect size. Data size is %s while indices size is %s"%(d,i)

        delegate(self, **kw)


    def density(self, *args, **kwargs):
        r"""Compute the density matrix of the Greens function

        Parameters
        ----------

        beta : float, optional
            Used for finite temperature density calculation with ``MeshReFreq``.

        Returns
        -------

        density_matrix : ndarray
            Single particle density matrix with shape ``target_shape``.

        Notes
        -----

        Only works for single mesh Greens functions with a, Matsubara,
        real-frequency, or Legendre mesh.
        """
        # Must be implemeted in a subclass
        raise NotImplementedError


    @property
    def rank(self):
        r"""int : The mesh rank (number of meshes)."""
        return self._rank

    @property
    def target_rank(self): 
        """int : The rank of the target space."""
        return self._target_rank

    @property
    def target_shape(self): 
        """(int, ...) : The shape of the target space."""
        return self._target_shape

    @property
    def mesh(self):
        """gf_mesh : The mesh of the Greens function."""
        return self._mesh

    @property
    def data(self):
        """ndarray : Raw data of the Greens function.

           Storage convention is ``self.data[x,y,z, ..., n0,n1,n2]``
           where ``x,y,z`` correspond to the mesh variables (the mesh) and 
           ``n0, n1, n2`` to the ``target_space``.
        """
        return self._data

    @property
    def indices(self):
        """GfIndices : The index object of the target space."""
        return self._indices

    def copy(self) : 
        """Deep copy of the Greens function.

        Returns
        -------
        G : Gf
            Copy of self.
        """
        return Gf (mesh = self._mesh.copy(), 
                   data = self._data.copy(), 
                   indices = self._indices.copy(), 
                   name = self.name)

    def copy_from(self, another):
        """Copy the data of another Greens function into self."""
        self._mesh.copy_from(another.mesh)
        assert self._data.shape == another._data.shape, "Shapes are incompatible: " + str(self._data.shape) + " vs " + str(another._data.shape)
        self._data[:] = another._data[:]
        self._indices = another._indices.copy()
        self.__check_invariants()

    def __repr__(self):
        return "Greens Function %s with mesh %s and target_rank %s: \n"%(self.name, self.mesh, self.target_rank)
 
    def __str__ (self): 
        return self.name if self.name else repr(self)

    #def __getitem__(self, key):
        #return self._data[key]

    #def __setitem__(self, key, val):
        #self[key] << val

    # -------------- Various operations -------------------------------------
    
    @property
    def real(self): 
        """Gf : A Greens function with a view of the real part."""
        return Gf(mesh = self._mesh, data = self._data.real, name = ("Re " + self.name) if self.name else '') 

    @property
    def imag(self): 
        """Gf : A Greens function with a view of the imaginary part."""
        return Gf(mesh = self._mesh, data = self._data.imag, name = ("Im " + self.name) if self.name else '') 
 
    def __lshift__(self, A):
        """ A can be two things:
          * G << any_init will init the GFBloc with the initializer
          * G << g2 where g2 is a GFBloc will copy g2 into self
        """
        if isinstance(A, Gf):
            if self is not A: # otherwise it is useless AND does not work !!
                assert self.mesh == A.mesh, "Green function meshes are not compatible:\n  %s\nand\n  %s" % (self.mesh, A.mesh)
                self.copy_from(A)
        return self
 
    def __le__(self, other): 
        raise RuntimeError(" Operator <= not defined ")

    # ---------- Addition   

    def __iadd__(self,arg):
        if isinstance(arg, Gf):
            assert type(self.mesh) == type(arg.mesh), "Can not add two Gf with meshes of different type"
            assert self.mesh == arg.mesh, "Can not add two Gf with different mesh"
            self._data += arg._data 
        else:
            raise NotImplemented
        return self

    def __add__(self,y):
        c = self.copy()
        c += y
        return c

    def __radd__(self,y): return self.__add__(y)

    # ---------- Substraction

    def __isub__(self,arg):
       if isinstance(arg, Gf):
           assert type(self.mesh) == type(arg.mesh), "Can not subtract two Gf with meshes of different type"
           assert self.mesh == arg.mesh, "Can not subtract two Gf with different mesh"
           self._data -= arg._data 
       else:
           raise NotImplemented
       return self

    def __sub__(self,y):
        c = self.copy()
        c -= y
        return c

    def __rsub__(self,y):
        c = (-1)*self.copy()
        c += y
        return c

   #----------------------------- other operations -----------------------------------

    def invert(self):
        """Inverts the Greens function (in place)."""

        if self.target_rank == 0: # Scalar target space
            self.data[:] = 1. / self.data
        elif self.target_rank == 2: # Matrix target space
            self.data[:] = np.linalg.inv(self.data)
        else:
            raise TypeError(
                "Inversion only makes sense for matrix or scalar_valued Greens functions")

    def inverse(self): 
        """Computes the inverse of the Greens function.

        Returns
        -------
        G : Gf (copy)
            The matrix/scalar inverse of the Greens function.
        """
        r = self.copy()
        r.invert()
        return r

    def transpose(self): 
        """Take the transpose of a matrix valued Greens function.

        Returns
        -------

        G : Gf (copy)
            The transpose of the Greens function.

        Notes
        -----

        Only implemented for single mesh matrix valued Greens functions.

        """

        # FIXME Why this assert ?
        assert self.rank == 1, "Transpose only implemented for single mesh Greens functions"
        assert self.target_rank == 2, "Transpose only implemented for matrix valued Greens functions"

        d = np.transpose(self.data.copy(), (0, 2, 1))
        return Gf(mesh = self.mesh, data= d, indices = self.indices.transpose())

    def conjugate(self):
        """Conjugate of the Greens function.

        Returns
        -------
        G : Gf (copy)
            Conjugate of the Greens function.
        """
        return Gf(mesh = self.mesh, data= np.conj(self.data), indices = self.indices)

    def zero(self):
        """Set all values to zero."""
        self._data[:] = 0

    def from_L_G_R(self, L, G, R):
        r"""Matrix transform of the target space of a matrix valued Greens function.

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

        Notes
        -----

        Only implemented for Greens functions with a single mesh.
        """

        assert self.rank == 1, "Only implemented for Greens functions with one mesh"
        assert self.target_rank == 2, "Matrix transform only valid for matrix valued Greens functions"

        assert len(L.shape) == 2, "L needs to be two dimensional"
        assert len(R.shape) == 2, "R needs to be two dimensional"

        assert L.shape[1] == G.target_shape[0], "Dimension mismatch between L and G"
        assert R.shape[0] == G.target_shape[1], "Dimension mismatch between G and R"

        assert L.shape[0] == self.target_shape[0], "Dimension mismatch between L and self"
        assert R.shape[1] == self.target_shape[1], "Dimension mismatch between R and self"

        self.data[...] = np.einsum('Mac,cd,db->Mab', L, G.data, R, optimize=True)

    def total_density(self, *args, **kwargs):
        """Compute total density.

        Returns
        -------
        density : float
            Total density of the Greens function.

        Notes
        -----
        Only implemented for single mesh Greens function with a,
        Matsubara, real-frequency, or Legendre mesh.

        """
        raise NotImplementedError

    #-----------------------------  IO  -----------------------------------
    
    def __reduce__(self):
        return call_factory_from_dict, (Gf, self.name, self.__reduce_to_dict__())

    def __reduce_to_dict__(self):
        d = {'mesh' : self._mesh, 'data' : self._data}
        if self.indices : d['indices'] = self.indices 
        return d

    _hdf5_format_ = 'Gf'

    @classmethod
    def __factory_from_dict__(cls, name, d):
        # Backward compatibility layer
        # Drop singularity from the element and ignore it
        d.pop('singularity', None)
        #
        r = cls(name = name, **d)
        # Backward compatibility layer
        # In the case of an ImFreq function, old archives did store only the >0
        # frequencies, we need to duplicate it for negative freq.
        # Same code as in the C++ h5_read for gf.
        need_unfold = isinstance(r.mesh, meshes.MeshImFreq) and r.mesh.positive_only() 
        return r if not need_unfold else wrapped_aux._make_gf_from_real_gf(r)
    
#---------------------------------------------------------

from h5.formats import register_class, register_backward_compatibility_method
register_class (Gf)

# A backward compatility function
def bckwd(hdf_scheme):
    # we know scheme is of the form GfM1_x_M2_s/tv3
    m, t= hdf_scheme[2:], '' # get rid of Gf
    for suffix in ['_s', 'Tv3', 'Tv4'] : 
        if m.endswith(suffix) :
            m, t = m[:-len(suffix)], suffix
            break
    return { 'mesh': 'Mesh'+m, 'indices': 'GfIndices'}

register_backward_compatibility_method("Gf", "Gf", bckwd)

