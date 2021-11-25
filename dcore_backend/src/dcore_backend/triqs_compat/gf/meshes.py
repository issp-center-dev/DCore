import numpy as np
from ..h5.archive import register_class

class Mesh(object):
    pass

    def __eq__(self, other):
        return type(self) == type(other) and self.hash == other.hash
    
    @property
    def points(self):
        return self._points

    @classmethod
    def n_points_fact(cls):
        """
        points.size/n_points
        """
        return 1
    
    @classmethod
    def default_n_points(cls):
        return -1

    @property
    def beta(self):
        return self._beta
    


class MeshReFreq(Mesh):
    """ Real frequency mesh """
    def __init__(self, omega_min, omega_max, n_points, beta=None):
        """

        Args:
            omega_min (float): min value of frequency
            omega_max (float): max value of frequency
            n_points (int): Number of frequencies
            beta (float): inverse temperature
        """
        self._omega_min = omega_min
        self._omega_max = omega_max
        self._n_max = n_points
        self._points = np.linspace(omega_min, omega_max, n_points)
        self._beta = beta

    @property
    def size(self):
        return self._points.size

    def __iter__(self):
        yield from self._points
    
    def values(self):
        return self._points

    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        group.create_group(key)
        group[key]['max'] = self._omega_max
        group[key]['min'] = self._omega_min
        group[key]['size'] = self.size
        group[key].write_attr('Format', 'MeshReFreq')

    @classmethod
    def __factory_from_dict__(cls, key, dict) :
        return cls(dict['min'], dict['max'], dict['size'])

class MeshImFreq(Mesh):
    """ Imaginary frequency mesh """
    def __init__(self, beta, n_points=None, statistic=None, S=None, n_max=None):
        """

        Args:
            beta (float): inverse temperature
            statistic (str): 'Fermion' or 'Boson'
            n_points (int):
                Number of non-negative frequencies
            S (str, optional): 'Fermion' or 'Boson'
            n_max (int): Old name for n_points
        """
        self._beta = beta
        if S is not None:
            assert statistic is None
            statistic = S
        if n_max is not None:
            assert n_points is None
            n_points = n_max
        assert isinstance(statistic, str), type(statistic)
        self._statistic = {'F': 'Fermion', 'B': 'Boson'}[statistic[0]]
        shift = 1 if self._statistic[0] == 'F' else 0
        self._points = 2*np.arange(-n_points, n_points) + shift
        self._values = 1J * (np.pi/self.beta) * self._points

    
    @classmethod
    def n_points_fact(cls):
        return 2

    @classmethod
    def default_n_points(cls):
        return 1025

    @property
    def hash(self):
        return hash(self._beta) + hash(self._statistic) + hash(self._points.tobytes())
    
    @property
    def statistic(self):
        return self._statistic
    
    @property
    def points(self):
        return self._points

    def values(self):
        return self._values

    @property
    def size(self):
        return self._points.size

    def __iter__(self):
        yield from self._values

    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        group.create_group(key)
        group[key]._raw_write('positive_freq_only', np.intc(0))
        group[key]['size'] = self.size
        group[key].create_group('domain')
        group[key]['domain']['beta'] = self.beta
        group[key]['domain']['statistic'] = self.statistic[0]
        group[key].write_attr('Format', 'MeshImFreq')
    
    @classmethod
    def __factory_from_dict__(cls, key, dict) :
        assert dict['positive_freq_only'] == False
        return cls(
            beta=dict['domain']['beta'],
            statistic=dict['domain']['statistic'],
            n_points=dict['size']//2
        )

class MeshImTime(Mesh):
    """ Imaginary time mesh """
    def __init__(self, beta, statistic, n_points):
        """

        Args:
            beta (float): inverse temperature
            statistic (str): 'Fermion' or 'Boson'
            n_points (int):
                Number of equidistant mesh points
        """
        assert isinstance(statistic, str)
        self._beta = beta
        self._statistic = {'F': 'Fermion', 'B': 'Boson'}[statistic[0]]
        self._points = np.linspace(0, beta, n_points)
    

    @property
    def statistic(self):
        return self._statistic
    
    @property
    def points(self):
        return self._points

    @property
    def size(self):
        return self._points.size

    def __iter__(self):
        yield from self._points

    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        pass
    
    @classmethod
    def __factory_from_dict__(cls, key, dict) :
        pass

class MeshLegendre(Mesh):
    """ Legendre mesh """
    def __init__(self, n_points):
        """

        Args:
            n_points (int): Number of Legendre polys
        """
        self._nl = n_points
        self._points = np.arange(n_points)

    @property
    def size(self):
        return self._nl

    def __iter__(self):
        yield from self._points

class MeshIR(Mesh):
    """ IR mesh """
    def __init__(self, basis):
        """

        Args:
            basis: IR basis 
        """
        self._points = np.arange(basis.size)
        self._basis = basis
        self._beta = basis.beta

    @property
    def size(self):
        return self._points.size

    def __iter__(self):
        yield from self._points

    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        group.create_group(key)
        group[key]['size'] = self.size
        group[key]['basis'] = self.basis
        group[key].write_attr('Format', 'MeshIR')

all_meshes = [MeshImFreq, MeshReFreq, MeshLegendre, MeshIR]

register_class(MeshImFreq)
register_class(MeshReFreq)
register_class(MeshLegendre)
register_class(MeshIR)