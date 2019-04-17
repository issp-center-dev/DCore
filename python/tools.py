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

import sys
import os
import numpy
import shlex
import subprocess
from itertools import *

from pytriqs.utility.h5diff import compare, failures
from pytriqs.archive.hdf_archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.operators import *

"""
THIS MODULE  MUST NOT DEPEND ON MPI!
"""

def h5diff(f1, f2, key, precision=1.e-6):
    """

    Modified version of pytriqs.utility.h5diff.h5diff
    key is the path of the data set to be compared: e.g., "dmft_out/Sigma_iw"

    """
    keys = key.split("/")
    h1 = HDFArchive(f1,'r')
    h2 = HDFArchive(f2,'r')

    for k in keys:
        h1 = h1[k]
        h2 = h2[k]

    compare(key, h1, h2, 0, precision)
    if failures :
        print ('-'*50, file=sys.stderr )
        print ('-'*20 + '  FAILED  ' +  '-'*20, file=sys.stderr)
        print ('-'*50, file=sys.stderr)
        for x in failures:
            print (x, file=sys.stderr)
            print ('-'*50, file=sys.stderr)
        raise RuntimeError("FAILED")

def to_spin_full_U_matrix(u_matrix):
    """
    To spin full four-index U matrix

    :param u_matrix: U_{ijkl}, the dimesion of each axis is the number of orbitals
    """
    n_orb = u_matrix.shape[0]

    u_matrix_spin_full = numpy.zeros((2, n_orb, 2, n_orb, 2, n_orb, 2, n_orb), dtype=u_matrix.dtype)

    for i1, i2 in product(range(2), repeat=2):
        u_matrix_spin_full[i1,:,i2,:,i1,:,i2,:] = u_matrix

    return u_matrix_spin_full.reshape((2*n_orb, 2*n_orb, 2*n_orb, 2*n_orb))


def make_real(bGf):
    for name, g in bGf:
        if type(g) == GfImFreq:
            bGf[name] = g.make_real_in_tau()
        elif type(g) == GfImTime:
            g.data[:,:,:] = g.data.real
            g.tail.data[:,:] = g.tail.data.real
        else:
            raise RuntimeError("Unsupported type " + str(type(g)))


def make_block_gf(gf_class, gf_struct, beta, n_points):
    """
    Make a BlockGf object

    :param gf_class: GfImFreq, GfImTime or GfLegendre
    :param gf_struct: structure of Green's function
    :param beta: inverse temperature
    :param n_points: number of points
    :return: object of BlockGf
    """

    assert isinstance(gf_struct, dict)

    blocks = []
    block_names = []
    for name, indices in gf_struct.items():
        assert isinstance(name, str)
        block_names.append(name)
        blocks.append(gf_class(indices=indices, beta=beta, n_points=n_points, name=name))
    return BlockGf(name_list=block_names, block_list=blocks, make_copies=True)


def compute_diag_basis(G0_iw):
    H_loc0 = {}
    for name, block in G0_iw:
        H_loc0[name] = numpy.array(block.tail[2])
    rot = {}
    for sp in H_loc0.keys():
        eigval, rot[sp] = numpy.linalg.eigh(H_loc0[sp])
    return rot

def gf_block_names(use_spin_orbit):
    """
    :param use_spin_orbit: bool
    :return:
        ['ud'] for use_spin_orbit = True, ['up', 'down'] for use_spin_orbit = False.
        This block structure is used for Green's functions throughout DCore.
    """

    if use_spin_orbit:
        return ['ud']
    else:
        return ['up', 'down']

def raise_if_mpi_imported():
    if 'pytriqs.utility.mpi' in sys.modules:
        raise RuntimeError("Error: MPI must not be imported in a non-MPI module! This indicates a bug in DCore.")

def convert_to_built_in_scalar_type(data):
    if isinstance(data, numpy.ndarray):
        return data
    elif numpy.issubdtype(type(data), numpy.integer):
        return int(data)
    elif numpy.issubdtype(type(data), numpy.float):
        return float(data)
    elif numpy.issubdtype(type(data), numpy.bool_):
        return bool(data)
    elif numpy.issubdtype(type(data), numpy.complex):
        return complex(data)


def symmetrize_spin(G):
    """
    Enforce time-reversal symmetry
    """
    # get spin labels, e.g., ['up', 'down']
    bnames = list(G.indices)
    assert len(bnames) == 2

    # average over spins
    G_ave = (G[bnames[0]] + G[bnames[1]])/2.
    G[bnames[0]] = G_ave.copy()
    G[bnames[1]] = G_ave.copy()

def launch_mpi_subprocesses(mpirun_command, rest_commands, output_file):
    """

    :param mpirun_command:  str
        e.g. "mpirun -np 8"
    :param rest_commands: list
        rest part of command

    """
    commands = shlex.split(mpirun_command)
    commands.extend(rest_commands)
    return_code = subprocess.call(commands, stdout=output_file, stderr=output_file)
    output_file.flush()
    if return_code:
        print("Error occurred while executing MPI program!")
        print("Return code: ", return_code)
        print("Command: ", ' '.join(commands))
        raise RuntimeError("Error occurred while executing MPI program! Output messages may be found in {}!".format(os.path.abspath(output_file.name)))

def extract_H0(G0_iw, hermitianize=True):
    """
    Extract non-interacting Hamiltonian elements from G0_iw
    """

    H0 = [numpy.array(block.tail[2]) for name, block in G0_iw]

    n_spin_orb = numpy.sum([b.shape[0] for b in H0])

    if G0_iw.n_blocks > 2:
        raise RuntimeError("n_blocks must be 1 or 2.")

    names = [name for name, block in G0_iw]

    data = numpy.zeros((n_spin_orb, n_spin_orb), dtype=complex)
    offset = 0
    for block in H0:
        block_dim = block.shape[0]
        data[offset:offset + block_dim, offset:offset + block_dim] = block
        offset += block_dim

    if hermitianize:
        data = 0.5 * (data.transpose().conj() + data)

    return data

# TODO: clean comments
def umat2dd(dcore_U):

    n_orb = dcore_U.shape[0]/2 # spin-1/2
    # n_orb = dcore_U.shape[0]
    print("n_orb:",n_orb)
    # extract density-density part
    dcore_U_len = len(dcore_U)
    alps_U = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_Uprime = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_J = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)

    # m_range = range(size)
    print("dcore_U_len:",dcore_U_len)
    for i, j in product(range(dcore_U_len), range(dcore_U_len)):
        alps_U[i, j] = dcore_U[i, j, i, j].real - dcore_U[i, j, j, i].real
        alps_Uprime[i, j] = dcore_U[i, j, i, j].real
        alps_J[i, j] = dcore_U[i, j, j, i].real
    #return alps_U, alps_Uprime, alps_J
    #def write_Umatrix(U, Uprime, J, norb):

    #Uout = numpy.zeros((n_orb, 2, n_orb, 2))
    Uout = numpy.zeros((2, n_orb, 2, n_orb))

    # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
    # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
    def func(u):
        print("shape:",u.shape)
        uout = u.reshape((2, n_orb, 2, n_orb))
        #uout = u.reshape((2, norb, 2, norb)).transpose(1, 0, 3, 2)
        return uout

    print("Uout:", Uout.shape, "alpsU:", alps_U.shape)
    U_four = func(alps_U)
    Uprime_four = func(alps_Uprime)
    J_four = func(alps_J)

    for a1, a2 in product(range(n_orb), repeat=2):
        for s1, s2 in product(range(2), repeat=2):  # spin-1/2
            if a1 == a2:
                Uout[s1, a1, s2, a2] = U_four[s1, a1, s2, a2]
                #Uout[a1, s1, a2, s2] = U_four[a1, s1, a2, s2]
            else:
                Uout[s1, a1, s2, a2] = Uprime_four[s1, a1, s2, a2] - J_four[s1, a1, s2, a2]
                #Uout[a1, s1, a2, s2] = Uprime_four[a1, s1, a2, s2] - J_four[a1, s1, a2, s2]

    Uout = Uout.reshape((2*n_orb, 2*n_orb))
    # with open('./Umatrix', 'w') as f:
    #    for i in range(2*norb):
    #        for j in range(2*norb):
    #            print('{:.15e} '.format(Uout[i, j].real), file=f, end="")
    #        print("", file=f)

# TODO: clean up
def dcore2alpscore(dcore_U):

    dcore_U_len = len(dcore_U)
    alps_U = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_Uprime = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_J = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)

    # m_range = range(size)
    for i, j in product(range(dcore_U_len), range(dcore_U_len)):
        alps_U[i, j] = dcore_U[i, j, i, j].real - dcore_U[i, j, j, i].real
        alps_Uprime[i, j] = dcore_U[i, j, i, j].real
        alps_J[i, j] = dcore_U[i, j, j, i].real
    return alps_U, alps_Uprime, alps_J

# TODO: clean up
def write_Umatrix(U, Uprime, J, norb):
    Uout = numpy.zeros((norb, 2, norb, 2))

    # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
    # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
    def func(u):
        uout = u.reshape((2, norb, 2, norb)).transpose(1, 0, 3, 2)
        return uout

    U_four = func(U)
    Uprime_four = func(Uprime)
    J_four = func(J)

    for a1, a2 in product(range(norb), repeat=2):
        for s1, s2 in product(range(2), repeat=2):  # spin-1/2
            if a1 == a2:
                Uout[a1, s1, a2, s2] = U_four[a1, s1, a2, s2]
            else:
                Uout[a1, s1, a2, s2] = Uprime_four[a1, s1, a2, s2] - J_four[a1, s1, a2, s2]

    Uout = Uout.reshape((2*norb, 2*norb))
    with open('./Umatrix', 'w') as f:
        for i in range(2*norb):
            for j in range(2*norb):
                print('{:.15e} '.format(Uout[i, j].real), file=f, end="")
            print("", file=f)
