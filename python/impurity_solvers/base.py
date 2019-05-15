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

from ..pytriqs_gf_compat import *
from pytriqs.archive import HDFArchive
from pytriqs.operators import *

from itertools import *
import numpy
import subprocess
import os
import shlex
import copy

from ..tools import *

class SolverBase(object):
    """
    This class define the common interface of all solvers.
    All solvers must inherit from this class!
    """

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : dict{str:list}
                    Structure of the Green's functions. It must be a
                    dictionary, which maps the name of each block of the
                    Green's function as a string to a list of integer
                    indices.
                    For example: ``{'up': [1,2,3], 'down', [1,2,3]}``.
                    Allowed block names are 'up', 'down' and 'ud'.
        u_mat numpy.ndarray:
            four-index spin-full U-matrix.
            The index of each axis runs as (up, orb1), (up, orb2), ..., (down, orb1), (down, orb2), ...
        n_iw : integer, optional
               Number of Matsubara frequencies used for the Green's functions.
               Used for input and output.
        """

        assert isinstance(gf_struct, dict)

        if len(gf_struct) > 2:
            raise RuntimeError("gf_struct cannot contain more than two blocks.")

        allowed_block_names = ['up', 'down', 'ud']
        if not numpy.all([name in allowed_block_names for name in gf_struct.keys()]):
            raise RuntimeError("Allowed block names are up, down, ud.")

        # Member variables
        self.beta = beta
        self.gf_struct = gf_struct
        self.u_mat = u_mat
        self.n_iw = n_iw
        self.dims = numpy.array([len(indices) for name, indices in self.gf_struct.items()])
        self.n_flavors = numpy.sum(self.dims)
        self.use_spin_orbit = (len(gf_struct) == 1)
        self.n_orb = int(self.n_flavors/2)

        self.block_names = gf_block_names(self.use_spin_orbit) # This should be consistent with gf_struct.

        if self.u_mat.shape[0] != 2*self.n_orb:
            raise RuntimeError("Size of U matrix is wrong!")

        # Input data
        # self.G0_iw must be set to proper values before solve() is called.
        self._G0_iw = make_block_gf(GfImFreq, gf_struct, beta, n_iw)

        # Output data.
        # These objects must be filled with results in solve().
        # Default value of Sigma_iw is 0, which will be used for the first iteration of DMFT self-consistent procedure.
        self._Sigma_iw = make_block_gf(GfImFreq, gf_struct, beta, n_iw)
        self._Sigma_iw.zero()
        self._Gimp_iw = make_block_gf(GfImFreq, gf_struct, beta, n_iw)
        self._Sigma_w = None
        #"self._G_l = None # Please define it if Legendre basis is used

    def name(self):
        return "Base solver"

    def set_G0_iw(self, new_G0_iw):
        self._G0_iw << new_G0_iw.copy()

    def get_Sigma_iw(self):
        return self._Sigma_iw.copy()

    def get_Gimp_iw(self):
        return self._Gimp_iw.copy()

    def get_Sigma_w(self):
        if self._Sigma_w is None:
            return None
        return self._Sigma_w.copy()

    def solve(self, rot, mpirun_command, params):
        """
        Solve the impudrity problem.

        self.G0_iw must be set before this function is called.

        INPUT:
          rot is a unitary matrix representing a basis transformation to reduce sign problem.
          If rot is None, no basis transformation is disabled.

          mpirun_command is a string (e.g. 'mpirun -np 8')

          params is a dict containing parameters.
          The following parameters must be set.
          Each impurity solver may require additional parameters.

              'random_seed_offset'  : int, off set value for random seed (mandatory)

              'calc_Sigma_w'        : bool, if real-frequency Sigma is computed (optional)
              'omega_min'           : float, lower limit of omega (must be set if calc_Sigma_iw is on).
              'omega_max'           : float, uppper limit of omega (must be set if calc_Sigma_iw is on).
              'n_omega'             : int, number of frequency points (must be set if calc_Sigma_iw is on).

        OUTPUT:
        ``Gimp_iw``, ``Sigma_iw`` must be calculated.

        """

        pass

        # Set self.Gimp_iw, self.G_tau, self.Sigma_iw

    @classmethod
    def is_gf_realomega_available(cls):
        return False


def creat_mapping_flatten_index(gf_struct):
    # Map (block_name, index) to an index in the flatten spin-orbital space
    # If blocks are 'up' and 'down', 'up' appears FIRST.
    to_flatten_index = {}
    if len(gf_struct) == 1:
        from_flatten_index = []
        for offset, index in enumerate(gf_struct['ud']):
            to_flatten_index[('ud', index)] = offset
            from_flatten_index.append(('ud', index))
    else:
        block_names = ['up', 'down']
        from_flatten_index = []
        offset = 0
        for name in block_names:
            for index in gf_struct[name]:
                to_flatten_index[(name, index)] = offset
                from_flatten_index.append((name, index))
                offset += 1

    return to_flatten_index, from_flatten_index


def make_h_int(u_mat, gf_struct):
    """
    Construct an operator representing the interacting Hamiltonian

    :param u_mat: four-index U matrix.
        The dimensions of each axis is spin * orbital.
           gf_struct: dict
    """

    n_orb = int(u_mat.shape[0]/2)
    _, from_flatten_index = creat_mapping_flatten_index(gf_struct)

    ham = Operator()
    for i1, i2, i3, i4 in product(range(2*n_orb), repeat=4):
        ham += 0.5 * u_mat[i1, i2, i3, i4] \
               * c_dag(*from_flatten_index[i1]) * c_dag(*from_flatten_index[i2]) \
               * c(*from_flatten_index[i4]) * c(*from_flatten_index[i3])

    return ham


def rotate_basis(rot, use_spin_orbit, u_matrix, Gfs=[], direction='forward'):
    """
    Rotate all Gf-like objects and U-matrix to the basis defined by rot

    :param direction:  'forward' or 'backward'
    :return u_matrix: rotated U matrix
    """

    if direction == 'forward':
        return _rotate_basis(rot, u_matrix, use_spin_orbit, Gfs)
    elif direction == 'backward':
        rot_conj_trans = {}
        for name, r in rot.items():
            rot_conj_trans[name] = r.conjugate().transpose()
        return _rotate_basis(rot_conj_trans, u_matrix, use_spin_orbit, Gfs)
    else:
        raise RuntimeError("Unknown direction " + direction)


def _rotate_basis(rot, u_matrix, use_spin_orbit, Gfs):
    """
    Rotate all Gf-like object and U matrix to a new local basis defined by "rot".
    :param rot: matrix
    :return u_matrix: rotated U matrix
    """

    if use_spin_orbit:
        rot_spin_full = rot['ud']
    else:
        n_orb = rot['up'].shape[0]
        rot_spin_full = numpy.zeros((2*n_orb, 2*n_orb), dtype=complex)
        rot_spin_full[0:n_orb, 0:n_orb] = rot['up']
        rot_spin_full[n_orb:, n_orb:] = rot['down']


    for G in Gfs:
        for bname, gf in G:
            gf.from_L_G_R(rot[bname].transpose().conjugate(), gf, rot[bname])

    if not u_matrix is None:
        return numpy.einsum("ijkl,im,jn,ko,lp", u_matrix,
                                    numpy.conj(rot_spin_full), numpy.conj(rot_spin_full), rot_spin_full, rot_spin_full)


class PytriqsMPISolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """

        Initialize the solver.

        """

        super(PytriqsMPISolver, self).__init__(beta, gf_struct, u_mat, n_iw)

    def _impl_module_name(self):
        return ""

    def solve(self, rot, mpirun_command, params_kw):

        params = copy.deepcopy(params_kw)

        # Write input parameters
        with HDFArchive(os.path.abspath('input.h5'), 'w') as h:
            h['beta'] = self.beta
            h['gf_struct'] = self.gf_struct
            h['u_mat'] = self.u_mat
            h['n_iw'] = self.n_iw
            if not rot is None:
                h['rot'] = rot
            h['G0_iw'] = self._G0_iw
            h['params'] = params

            if 'calc_Sigma_w' in params and params['calc_Sigma_w']:
                h['calc_Sigma_w'] = True
                for k in ['omega_min', 'omega_max', 'n_omega']:
                    h[k] = params[k]
            else:
                h['calc_Sigma_w'] = False

        # Run a working horse
        commands = [sys.executable, "-m", self._impl_module_name()]
        commands.append(os.path.abspath('./input.h5'))
        commands.append(os.path.abspath('./output.h5'))

        with open('./output', 'w') as output_file:
            launch_mpi_subprocesses(mpirun_command, commands, output_file)
        with open('./output', 'r') as output_file:
            for line in output_file:
                print(line)

        # Read results
        with HDFArchive(os.path.abspath('output.h5'), 'r') as h:
            self._Sigma_iw << h['Sigma_iw']
            self._Gimp_iw << h['Gimp_iw']
            if 'Sigma_w' in h:
                self._Sigma_w = h['Sigma_w']

    def name(self):
        return "PytriqsMPISolver"
