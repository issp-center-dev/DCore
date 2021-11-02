import numpy as np

class Mesh(object):
   pass

class MeshReFreq(Mesh):
    """ Real frequency mesh """
    def __init__(self, omega_min=None, omega_max=None, n_max=None):
        """

        Args:
            omega_min (float): min value of frequency
            omega_max (float): max value of frequency
            n_max (int): Number of frequencies
        """
        self._omega_min = omega_min
        self._omega_max = omega_max
        self._n_max = n_max
        self._points = np.linspace(omega_min, omega_max, n_max)

    @property
    def size(self):
        return self._points.size

    def __iter__(self):
        yield from self._points

class MeshImFreq(Mesh):
    """ Imaginary frequency mesh """
    def __init__(self, beta, *kwargs):
        """

        Args:
            beta (float): inverse temperature
            S (str): 'Fermion' or 'Boson'
            n_points (int):
                Number of non-negative frequencies
            points (ndarray):
                integers representation of freqquencies
        """
        self._beta = beta
        self._statistic = kwargs['S']
        assert isinstance(self._statistic, str)
        if 'n_points' in kwargs:
            self._points = np.arange(-kwargs['n_points'], kwargs['n_points'])
        else:
            raise RuntimeError("points was not properly set.")
    
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

all_meshes = [MeshImFreq, MeshReFreq]
