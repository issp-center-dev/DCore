from copy import deepcopy, copy
import numpy as np
from .meshes import MeshImFreq
from dcore.triqs_compat.h5.archive import vls_dt
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

    """
    def __init__(self, **kw): # enforce keyword only policy 
        def delegate(self, data=None, name='', beta=None, statistic='Fermion', mesh=None):
            self.data = data
            self.target_shape = self.data[1:]
            self.name = name
            self.beta = beta
            self.statistic = statistic
            if mesh is None:
                mesh = MeshImFreq(beta, statistic, self.data.shape[0]//2)
            self.mesh = mesh

        delegate(self, **kw)
    

    def __lshift__(self, A):
        """ Substitute a new gf object (copy) """
        if isinstance(A, Gf):
            for name in ['data', 'target_shape', 'name', 'beta', 'statistic']:
                self.__setattr__(name, copy(A.__getattribute__(name)))
        elif isinstance(A, np.ndarray):
            self.data[...] = A
        else:
            raise RuntimeError("Invalid type of A!")

    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        group.create_group(key)
        group[key].write_attr('Format', 'Gf')
        group[key]['data'] = self.data
        group[key]['mesh'] = self.mesh
        assert self.data.ndim == 3
        #group[key].create_group('indices')
        #group[key]['indices'].write_attr('Format', 'GfIndices')
        #group[key]['indices']['left'] = np.array([str(i) for i in range(self.data.shape[1])], dtype=vls_dt)
        #group[key]['indices']['right'] = np.array([str(i) for i in range(self.data.shape[2])], dtype=vls_dt)
        ## Write Indices
        #__h5_writers[type(self.mesh)](self, group, key)
        #if isinstance(self.mesh, MeshImFreq):
            #__write_hdf5__as_gfimfreq(g, group, key):
        #else:
            #raise RuntimeError(f"Mesh {type(self.mesh)} is not writable!")

#def __write_hdf5__as_gfimfreq(g, group, key):
    #pass

#class GfImFreq(object): 
    #"""
    #Imaginary-frequency Green's function object
    #"""
    #def __init__(self, **kw): # enforce keyword only policy 
        #super().__init__(**kw)
    #
    #def __write_hdf5__(self, group, key):
        #""" Write to a HDF5 file"""
        #pass
#
    #@property
    #def n_points(self):
        #return self.data.shape[0]//2

#__h5_writers = {
    #MeshImFreq:  __write_hdf5__as_gfimfreq
#}