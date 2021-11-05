import numpy as np
from ..h5.archive import register_class

class Mesh(object):
   pass

class MeshReFreq(Mesh):
    """ Real frequency mesh """
    def __init__(self, omega_min, omega_max, n_points):
        """

        Args:
            omega_min (float): min value of frequency
            omega_max (float): max value of frequency
            n_points (int): Number of frequencies
        """
        self._omega_min = omega_min
        self._omega_max = omega_max
        self._n_max = n_points
        self._points = np.linspace(omega_min, omega_max, n_points)

    @property
    def size(self):
        return self._points.size

    def __iter__(self):
        yield from self._points

class MeshImFreq(Mesh):
    """ Imaginary frequency mesh """
    def __init__(self, beta, statistic, n_points):
        """

        Args:
            beta (float): inverse temperature
            statistic (str): 'Fermion' or 'Boson'
            n_points (int):
                Number of non-negative frequencies
        """
        assert isinstance(statistic, str)
        self._beta = beta
        self._statistic = {'F': 'Fermion', 'B': 'Boson'}[statistic[0]]
        self._points = np.arange(-n_points, n_points)
    
    @property
    def beta(self):
        return self._beta
    
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
        group.create_group(key)
        group[key]['positive_freq_only'] = False
        group[key]['size'] = self.size
        group[key].create_group('domain')
        group[key]['domain']['beta'] = self.beta
        group[key]['domain']['statistic'] = self.statistic[0]
        group[key].write_attr('Format', 'MeshImFreq')
    
    @classmethod
    def __factory_from_dict__(cls, key, dict) :
        assert dict['positive_freq_only'] == False
        return cls(
            dict['domain']['beta'],
            dict['domain']['statistic'],
            dict['size']//2
        )

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

all_meshes = [MeshImFreq, MeshReFreq, MeshLegendre]

register_class(MeshImFreq)
register_class(MeshReFreq)
register_class(MeshLegendre)
