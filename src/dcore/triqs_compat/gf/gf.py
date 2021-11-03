from copy import deepcopy, copy
import numpy as np
from .meshes import MeshImFreq
from dcore.triqs_compat.h5.archive import vls_dt, register_class
import h5py

def _to_fixed_length_utf8_array(str_list):
    #str_list = [(s.decode(encoding='utf-8') if isinstance(s, bytes) else s) for s in str_list]
    #print("debug", str_list)
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
            left_indices = list(map(str, np.arange(self.data.shape[1])))
            right_indices = list(map(str, np.arange(self.data.shape[2])))
            self.indices = GfIndices([left_indices, right_indices])

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

register_class(GfIndices)
register_class(Gf)