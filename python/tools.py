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
import ast

from pytriqs.utility.h5diff import compare, failures
from pytriqs.utility.h5diff import h5diff as h5diff_org
from pytriqs.archive.hdf_archive import HDFArchive
from .pytriqs_gf_compat import *
from pytriqs.operators import *
import scipy

from pytriqs import version

triqs_major_version = int(version.version.split('.')[0])

"""
THIS MODULE MUST NOT DEPEND ON MPI!
"""

def h5diff(f1, f2, key=None, precision=1.e-6):
    """

    Modified version of pytriqs.utility.h5diff.h5diff
    key is the path of the data set to be compared: e.g., "dmft_out/Sigma_iw"

    """
    if key is None:
        h5diff_org(os.path.abspath(f1), os.path.abspath(f2), precision)
        return

    keys = key.split("/")
    h1 = HDFArchive(os.path.abspath(f1),'r')
    h2 = HDFArchive(os.path.abspath(f2),'r')

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
        indices_str = list(map(str, indices))
        blocks.append(gf_class(indices=indices_str, beta=beta, n_points=n_points, name=name))
    return BlockGf(name_list=block_names, block_list=blocks, make_copies=True)

def get_block_size(g):
    """
    Get block size of g

    Parameters
    ----------
    g : object of GfImFreq, GfImTime or GfLegendre

    Returns
    -------
    block_size : int
        block of g (e.g. number of orbitals)
    """

    if triqs_major_version == 1:
        assert isinstance(g, GfImFreq) or isinstance(g, GfImTime) or isinstance(g, GfLegendre), 'Unsupported type {}'.format(type(g))
        return len(g.indices)
    elif triqs_major_version >= 2:
        assert len(g.indices[0]) == len(g.indices[1])
        return len(g.indices[0])

def extract_H0_from_tail(G0_iw):
    if isinstance(G0_iw, BlockGf):
        return {name:extract_H0_from_tail(b) for name, b in G0_iw}
    elif isinstance(G0_iw.mesh, MeshImFreq):
        if triqs_major_version == 1:
            return numpy.array(G0_iw.tail[2])
        elif triqs_major_version >= 2:
           import pytriqs.gf.gf_fnt as gf_fnt
           assert len(G0_iw.target_shape) in [0,2], "extract_H0_from_tail(G0_iw) requires a matrix or scalar_valued Green function"
           assert gf_fnt.is_gf_hermitian(G0_iw), "extract_H0_from_tail(G0_iw) requires a Green function with the property G0_iw[iw][i,j] = conj(G0_iw[-iw][j,i])"
           tail, err = gf_fnt.fit_hermitian_tail(G0_iw)
           if err > 1e-5:
               print("WARNING: delta extraction encountered a sizeable tail-fit error: ", err)
           return tail[2]
    else:
        raise RuntimeError('extract_H0_from_tail does not support type {}'.format(type(G0_iw)))

def compute_diag_basis(G0_iw):
    H_loc0 = extract_H0_from_tail(G0_iw)
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

def extract_H0(G0_iw, block_names, hermitianize=True):
    """
    Extract non-interacting Hamiltonian elements from G0_iw
    """

    assert isinstance(block_names, list)

    H0_dict  = extract_H0_from_tail(G0_iw)
    H0 = [H0_dict[b] for b in block_names]

    n_spin_orb = numpy.sum([b.shape[0] for b in H0])

    if G0_iw.n_blocks > 2:
        raise RuntimeError("n_blocks must be 1 or 2.")

    data = numpy.zeros((n_spin_orb, n_spin_orb), dtype=complex)
    offset = 0
    for block in H0:
        block_dim = block.shape[0]
        data[offset:offset + block_dim, offset:offset + block_dim] = block
        offset += block_dim

    if hermitianize:
        data = 0.5 * (data.transpose().conj() + data)

    return data


def umat2dd(dcore_U):

    n_orb = dcore_U.shape[0]/2  # spin-1/2

    # extract density-density part
    dcore_U_len = len(dcore_U)
    Uout = numpy.zeros((2*n_orb, 2*n_orb, 2*n_orb, 2*n_orb))

    for i, j, k, l in product(range(dcore_U_len), range(dcore_U_len), range(dcore_U_len), range(dcore_U_len)):
        if (i == k and j == l) or (i == l and j == k):
            Uout[i, j, k, l] = dcore_U[i, j, k, l]

    return Uout

def pauli_matrix():
    pauli_mat = []
    pauli_mat.append(numpy.array([[0, 1], [1,0]], dtype=complex))
    pauli_mat.append(numpy.array([[0, -1J], [1J, 0]], dtype=complex))
    pauli_mat.append(numpy.array([[1, 0], [0, -1]], dtype=complex))
    return pauli_mat

def spin_moments_sh(dm_sh):
    """
    Compute spin moments on shells.
    dm_sh must contain density matrices.
    A fully polarized S=1/2 spin gives 1/2.
    """

    pauli_mat = pauli_matrix()

    assert numpy.allclose(numpy.dot(pauli_mat[0],pauli_mat[1]), 1J*pauli_mat[2])

    spin_moments = []
    for ish in range(len(dm_sh)):
        dm_dict = dm_sh[ish]
        if len(dm_dict) == 1:
            dm = dm_dict['ud']
        else:
            dm = scipy.linalg.block_diag(dm_dict['up'], dm_dict['down'])

        assert dm.shape[0] == dm.shape[1]

        norb = dm.shape[0]//2
 
        dm = dm.reshape((2, norb, 2, norb))

        s = numpy.array([0.5*numpy.einsum('st, sntn', pauli_mat[i], dm).real for i in range(3)])
        spin_moments.append(s)

    return spin_moments


def read_potential(filename, mat):
    if not os.path.exists(filename):
        print("Error: file '{}' not found".format(filename))
        exit(1)
    print("Reading '{}'...".format(filename))

    filled = numpy.full(mat.shape, False)
    try:
        with open(filename, 'r') as f:
            for line in f:
                # skip comment line
                if line[0] == '#':
                    continue

                array = line.split()
                assert len(array) == 5
                sp = int(array[0])
                o1 = int(array[1])
                o2 = int(array[2])
                val = complex(float(array[3]), float(array[4]))
                if filled[sp, o1, o2]:
                    raise Exception("duplicate components: " + line)
                mat[sp, o1, o2] = val
                filled[sp, o1, o2] = True
    except Exception as e:
        print("Error:", e)
        print(line, end="")
        exit(1)


def set_potential(input_str, name, n_inequiv_shells, dim_sh, spin_orbit):
    """

    Parameters
    ----------
    input_str
    name
    n_inequiv_shells
    dim_sh
    spin_orbit

    Returns
    -------
    numpy.ndarray with complex type

        shape = (2, norb, norb)     w/  spin-orbit
                (1, 2*norb, 2*norb) w/o spin-orbit

    """

    print("\nInterpreting {} = {}".format(name, repr(input_str)))

    # init potential matrix
    if spin_orbit:
        pot = [numpy.zeros((1, dim_sh[ish], dim_sh[ish]), numpy.complex_) for ish in range(n_inequiv_shells)]
    else:
        pot = [numpy.zeros((2, dim_sh[ish], dim_sh[ish]), numpy.complex_) for ish in range(n_inequiv_shells)]

    # read potential matrix
    if input_str != 'None':
        try:
            files = ast.literal_eval(input_str)
            assert isinstance(files, dict), "should be dictionary"
            assert all([ish < n_inequiv_shells for ish in files.keys()]), "The keys must fulfill: key < n_inequiv_shells"
        except Exception as e:
            print("Error: %s =" % name, input_str)
            print(e)
            exit(1)

        for ish, file in files.items():
            read_potential(file, pot[ish])

    # print potential
    print("\n--- results for %s" % name)
    for ish, pot_ish in enumerate(pot):
        print("ish =", ish)
        for sp in range(pot_ish.shape[0]):
            print("sp =", sp)
            print(pot_ish[sp])

    return pot


def save_Sigma_iw_sh_txt(filename, Sigma_iw_sh, spin_names):
    """
    Save a list of self-energy in a text file

    :param filename:
    :param Sigma_iw_sh:
    :param spin_names: ['up', 'down'] or ['ud']
    """

    n_sh = len(Sigma_iw_sh)

    assert isinstance(spin_names, list)

    with open(filename, 'w') as fo:
        print("# Local self energy at imaginary frequency", file=fo)
        #
        # Column information
        #
        print("# [Column] Data", file=fo)
        print("# [1] Frequency", file=fo)
        icol = 1
        for ish in range(n_sh):
            for sp in spin_names:
                block_dim = Sigma_iw_sh[ish][sp].data.shape[1]
                for iorb, jorb in product(range(block_dim), repeat=2):
                    icol += 1
                    print("# [%d] Re(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, sp, iorb, jorb), file=fo)
                    icol += 1
                    print("# [%d] Im(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, sp, iorb, jorb), file=fo)

        #
        # Write data
        #
        omega = [x for x in Sigma_iw_sh[0].mesh]
        for iom in range(len(omega)):
            print("%f " % omega[iom].imag, end="", file=fo)
            for ish in range(n_sh):
                for sp in spin_names:
                    block_dim = Sigma_iw_sh[ish][sp].data.shape[1]
                    for iorb, jorb in product(range(block_dim), repeat=2):
                        print("%f %f " % (Sigma_iw_sh[ish][sp].data[iom, iorb, jorb].real,
                                          Sigma_iw_sh[ish][sp].data[iom, iorb, jorb].imag), end="", file=fo)
            print("", file=fo)


def readline_ignoring_comment(f):
    """
    Read a line ignoring lines starting with '#' or '%'

    :param f:
    :return:
    """

    while True:
        line = f.readline()
        is_comment = line.startswith('#')
        if not line or not is_comment:
            break
    return line


def load_Sigma_iw_sh_txt(filename, Sigma_iw_sh, spin_names, atol_omega=1e-2):
    """
    Load a list of self-energy from a text file.
    Tails will be set to zero.

    :param filename:
    :param Sigma_iw_sh: All elements must be allocated to correct shapes
    :param spin_names: ['up', 'down'] or ['ud']
    :param atol: absolute torelance for imaginary frequencies
    """

    n_sh = len(Sigma_iw_sh)

    assert isinstance(spin_names, list)

    data = numpy.loadtxt(filename)

    omega_imag = data[:, 0]
    nomega = len(omega_imag)
    if not numpy.allclose(omega_imag, numpy.array([complex(x) for x in Sigma_iw_sh[0].mesh]).imag, rtol=1, atol=atol_omega):
        raise RuntimeError("Mesh is not compatible!")

    for iom in range(nomega):
        icol = 1
        for ish in range(n_sh):
            for sp in spin_names:
                # FIXME: How to set zero to tail in TRIQS 2.x?
                #Sigma_iw_sh[ish][sp].tail.zero()
                block_dim = Sigma_iw_sh[ish][sp].data.shape[1]
                for iorb, jorb in product(range(block_dim), repeat=2):
                    re = data[iom, icol]
                    icol += 1
                    imag = data[iom, icol]
                    icol += 1
                    Sigma_iw_sh[ish][sp].data[iom, iorb, jorb] = complex(re, imag)

def __complex_to_float_array(a):
    return a.view(float).reshape(a.shape + (2,))

def __float_to_complex_array(a):
    return a.view(complex).reshape(a.shape[:-1])

def save_giw(h5file, path, g):
    """

    Save an object of GfImFreq to a HDF5 file

    Parameters
    ----------
    h5file : HDF5 archive (h5py)
    path : str
    g : GfImFreq
      object to be saved

    """

    assert isinstance(g, GfImFreq), 'Type {} is not supported by save_giw'.format(type(g))

    h5file[path + '/__version'] = 'DCore_GfImFreq_v1'
    h5file[path + '/data'] = __complex_to_float_array(g.data)
    h5file[path + '/tail'] = __complex_to_float_array(g.tail.data)
    h5file[path + '/wn'] = numpy.array([x for x in g.mesh]).imag


def load_giw(h5file, path, g):
    """

    Load data to an object of GfImFreq from a HDF5 file
    The mesh and shape of g must be set in advance.

    Parameters
    ----------
    h5file : HDF5 archive (h5py)
    path : str
    g : GfImFreq
      data will be loaded to g

    """

    assert isinstance(g, GfImFreq)
    assert h5file[path + '/__version'][()] == 'DCore_GfImFreq_v1'

    g.data[...] = __float_to_complex_array(h5file[path + '/data'][()])
    g.tail.data[...] = __float_to_complex_array(h5file[path + '/tail'][()])

    omega_imag = numpy.array([complex(x) for x in g.mesh]).imag
    if not numpy.allclose(omega_imag, h5file[path + '/wn'][()]):
        raise RuntimeError("Mesh is not compatible!")


def make_empty_dir(dir_path):
    """

    Prepare a working directory as an empty directory.
    Any existing file or directory at the specified path is removed.

    Parameters
    ----------
    dir_path : str
        Path to a directory

    """
    import shutil

    if os.path.exists(dir_path):
        if os.path.isfile(dir_path):
            os.remove(dir_path)
        elif os.path.isdir(dir_path):
            shutil.rmtree(dir_path)
        else:
            raise RuntimeError("Failed to remove " + dir_path)

    os.makedirs(dir_path)


def make_hermite_conjugate(Sigma_iw, check_only=False):
    """
    Make Sigma(iw_n) or G(iwn_n) hermite
    Return max difference.
    """
    flag = True
    max_diff = 0.0
    for name, g in Sigma_iw:
        # symmetrize tail
        for i in range(g.tail.data.shape[0]):
            g.tail.data[i, :, :] = 0.5 * (g.tail.data[i, :, :] + g.tail.data[i, :, :].conjugate().transpose())

        n_points = g.data.shape[0]//2
        for i in range(n_points):
            diff = numpy.amax(numpy.abs(g.data[i + n_points, :, :]-g.data[n_points - i - 1, :, :].conj().transpose()))
            max_diff = max(max_diff, diff)
            if not check_only:
                g.data[i + n_points, :, :] = g.data[n_points - i - 1, :, :].conj().transpose()
    return max_diff


