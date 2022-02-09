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


import sys
import os
import numpy
import shlex
import subprocess
from itertools import *
import ast
import math
import shutil
import scipy
from scipy import linalg as scipy_linalg

from dcore._dispatcher import h5diff as h5diff_org, compare, \
    BlockGf, HDFArchive, failures, MeshImFreq, fit_hermitian_tail, Gf, GfImFreq

"""
THIS MODULE MUST NOT DEPEND ON MPI!
"""


def h5diff(f1, f2, key=None, precision=1.e-6):
    """

    Modified version of triqs.utility.h5diff.h5diff
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

    for i1, i2 in product(list(range(2)), repeat=2):
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
    for name, indices in list(gf_struct.items()):
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

    assert len(g.indices[0]) == len(g.indices[1])
    return len(g.indices[0])

def extract_H0_from_tail(G0_iw):
    if isinstance(G0_iw, BlockGf):
        return {name:extract_H0_from_tail(b) for name, b in G0_iw}
    elif isinstance(G0_iw.mesh, MeshImFreq):
       assert len(G0_iw.target_shape) in [0,2], "extract_H0_from_tail(G0_iw) requires a matrix or scalar_valued Green function"
       tail, err = fit_hermitian_tail(G0_iw)
       if err > 1e-5:
           print("WARNING: delta extraction encountered a sizeable tail-fit error: ", err)
       return tail[2]
    else:
        raise RuntimeError('extract_H0_from_tail does not support type {}'.format(type(G0_iw)))

def compute_diag_basis(G0_iw):
    H_loc0 = extract_H0_from_tail(G0_iw)
    rot = {}
    for sp in list(H_loc0.keys()):
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
    if 'triqs.utility.mpi' in sys.modules or 'mpi4py' in sys.modules:
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
    raise_if_mpi_imported()
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

def _ph_symmetrize(eps):
    N = eps.size
    half_eps = eps[0:N//2]
    if N%2 == 0:
        return numpy.hstack((half_eps, -half_eps[::-1]))
    else:
        return numpy.hstack((half_eps, 0.0, -half_eps[::-1]))



def fit_delta_iw(delta_iw, beta, n_bath, n_fit, ph_symmetric, verbose, **fit_params):
    """
    Fit Delta(iw) using scipy

    scipy.optimize.fmin_bfgs
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_bfgs.html

    Parameters
    ----------
    delta_iw: [numpy.ndarray] (n_w, n_orb, n_orb)
    beta: [float] 1/T
    n_bath: [int] number of bath
    n_fit: [int] number of repetition of fitting
    ph_symmetric: [bool] particle-hole symmetric
    **fit_params: [dict] optional parameters to the fitting function

    Returns
    -------
    eps: [numpy.ndarray] (n_bath,) bath levels
    hyb: [numpy.ndarray] (n_orb, n_bath) hybridization parameters

    """
    from scipy import optimize

    n_w = delta_iw.shape[0]
    n_orb = delta_iw.shape[1]
    assert delta_iw.shape[2] == n_orb

    # fermionic Matsubara freqs
    freqs = numpy.array([1j * (2*i+1) * math.pi / beta for i in range(n_w)])

    # Define distance between delta_iw and delta_fit
    # delta_fit = sum_{l=1}^{n_bath} V_{o1, l} * V_{l, o2} / (iw - eps_{l})
    def distance(x):
        _eps = x[0:n_bath]
        if ph_symmetric:
            _eps = _ph_symmetrize(_eps)
        _hyb = x[n_bath:].reshape(n_orb, n_bath)

        # denom[i,j] = (freqs[i] - eps[j])
        denom = freqs[:, None] - _eps[None, :]

        # sum over bath index l
        delta_fit = numpy.einsum('al, bl, wl->wab', _hyb, _hyb.conj(), numpy.reciprocal(denom))

        # squared error
        return numpy.square(numpy.linalg.norm(delta_iw - delta_fit))

    # Determine eps and V which minimize the distance between delta_iw and delta_fit
    dis_min = 1.0e+10
    # [0:n_bath] -> eps_{l},  [n_bath:n_bath+n_orb*n_bath] -> V_{o,l}
    result_best = numpy.zeros(n_bath + n_orb * n_bath, dtype=float)
    for l in range(n_fit):
        # initial guess, random values in the range [-1:1]
        x0 = 2 * numpy.random.rand(result_best.size) - 1

        # fitting
        result = optimize.fmin_bfgs(distance, x0, **fit_params)
        if(verbose):
            print(" ", result)
        dis = distance(result)

        # update result_best
        if dis < dis_min:
            dis_min = dis
            result_best = result.copy()

    eps = result_best[0:n_bath]
    if ph_symmetric:
        eps = _ph_symmetrize(eps)
    hyb = result_best[n_bath:].reshape(n_orb, n_bath)
    return eps, hyb


def extract_bath_params(delta_iw, beta, block_names, n_bath, ph_symmetric=False, n_fit=5, fit_gtol=1e-5, verbose=False):
    """
    Determine bath parameters by fitting Delta(iw)

    Parameters
    ----------
    delta_iw: [block Gf] Delta(iw)
    beta: [float] 1/T
    block_names: [list] block names
    n_bath: [int] number of bath
    ph_symmetric: [bool] particle-hole symmetric
    n_fit: [int] number of repetition of fitting. The best fit result will be taken.
    fit_gtol: [float] A fitting parameter: Gradient norm must be less than gtol before successful termination.

    Returns
    -------
    eps: [numpy.ndarray] (2*n_bath,) bath levels
    hyb: [numpy.ndarray] (2*n_orb, 2*n_bath) hybridization parameters

    """

    n_orb = delta_iw[block_names[0]].data.shape[1]
    n_blocks = len(block_names)

    # These arrays will be returned
    eps_full = numpy.zeros((n_bath * n_blocks,), dtype=float)
    hyb_full = numpy.zeros((n_orb * n_blocks, n_bath * n_blocks), dtype=float)

    if n_bath == 0:
        return eps_full, hyb_full

    # fitting parameters
    fit_params = {
        "gtol": fit_gtol,
        "disp": verbose,
    }

    print("\nDetermine bath parameters by fitting Delta(iw)")
    for key, val in list(fit_params.items()):
        print("  {} : {}".format(key, val))

    # bath parameters for each block
    eps_list = []
    hyb_list = []
    for b in block_names:
        # fit Delta(iw)
        if(verbose):
            print("\nblock =", b)
        # data.shape == (n_w, n_orb, n_orb)
        n_w = delta_iw[b].data.shape[0]
        assert delta_iw[b].data.shape[1] == delta_iw[b].data.shape[2] == n_orb

        # use only positive Matsubara freq
        eps, hyb = fit_delta_iw(delta_iw[b].data[n_w//2:n_w, :, :], beta, n_bath, n_fit, ph_symmetric, verbose, **fit_params)
        assert eps.shape == (n_bath,)
        assert hyb.shape == (n_orb, n_bath)
        eps_list.append(eps)
        hyb_list.append(hyb)

    # Combine spin blocks
    # eps_full = {eps[up], eps[dn]}
    for i, block in enumerate(eps_list):
        n = block.shape[0]
        eps_full[n*i:n*(i+1)] = block

    # hyb_full = {{hyb[up], 0}, {0, hyb[dn]}}
    for i, block in enumerate(hyb_list):
        m, n = block.shape
        hyb_full[m*i:m*(i+1), n*i:n*(i+1)] = block

    print("\nfitting results")
    print("  eps[l]    hyb[0,l]  hyb[1,l]  ...")
    for l in range(eps_full.size):
        print(" %9.5f" %eps_full[l], end="")
        for orb in range(hyb_full.shape[0]):
            print(" %9.5f" %hyb_full[orb, l], end="")
        print("")

    return eps_full, hyb_full


def umat2dd(dcore_U):

    n_orb = dcore_U.shape[0]//2  # spin-1/2

    # extract density-density part
    dcore_U_len = len(dcore_U)
    Uout = numpy.zeros((2*n_orb, 2*n_orb, 2*n_orb, 2*n_orb), dtype=complex)

    for i, j, k, l in product(list(range(dcore_U_len)), repeat=4):
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
            dm = scipy_linalg.block_diag(dm_dict['up'], dm_dict['down'])

        assert dm.shape[0] == dm.shape[1]

        norb = dm.shape[0]//2

        dm = dm.reshape((2, norb, 2, norb))

        s = numpy.array([0.5*numpy.einsum('st, sntn', pauli_mat[i], dm).real for i in range(3)])
        spin_moments.append(s)

    return spin_moments


def read_potential(filename, mat):
    if not os.path.exists(filename):
        print("Error: file '{}' not found".format(filename), file=sys.stderr)
        exit(1)
    print("Reading '{}'...".format(filename))

    filled = numpy.full(mat.shape, False)
    try:
        with open(filename, 'r') as f:
            for line in f:
                line_comment_removed = line.split('#')[0]  # remove comment
                array = line_comment_removed.split()
                if len(array) == 0:  # skip an empty line
                    continue
                if len(array) != 5:
                    raise Exception(f"expect 5 columns, but {len(array)} columns entered")
                sp = int(array[0])
                o1 = int(array[1])
                o2 = int(array[2])
                val = complex(float(array[3]), float(array[4]))
                if filled[sp, o1, o2]:
                    raise Exception("duplicate components: " + line)
                mat[sp, o1, o2] = val
                filled[sp, o1, o2] = True
    except Exception as e:
        print("Error:", e, file=sys.stderr)
        print(line, end="", file=sys.stderr)
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
            assert all([ish < n_inequiv_shells for ish in list(files.keys())]), "The keys must fulfill: key < n_inequiv_shells"
        except Exception as e:
            print("Error: %s =" % name, input_str, file=sys.stderr)
            print(e, file=sys.stderr)
            exit(1)

        for ish, file in list(files.items()):
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
                for iorb, jorb in product(list(range(block_dim)), repeat=2):
                    icol += 1
                    print("# [%d] Re(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, sp, iorb, jorb), file=fo)
                    icol += 1
                    print("# [%d] Im(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, sp, iorb, jorb), file=fo)

        #
        # Write data
        #
        omega = [x for x in Sigma_iw_sh[0].mesh]
        for iom in range(len(omega)):
            print("%.15e " % omega[iom].imag, end="", file=fo)
            for ish in range(n_sh):
                for sp in spin_names:
                    block_dim = Sigma_iw_sh[ish][sp].data.shape[1]
                    for iorb, jorb in product(list(range(block_dim)), repeat=2):
                        print("%.15e %.15e " % (Sigma_iw_sh[ish][sp].data[iom, iorb, jorb].real,
                                          Sigma_iw_sh[ish][sp].data[iom, iorb, jorb].imag), end="", file=fo)
            print("", file=fo)


def save_Sigma_w_sh_txt(filename, Sigma_w_sh, spin_names):
    """
    Save a list of self-energy in a text file (real frequency ver.)

    :param filename:
    :param Sigma_w_sh:
    :param spin_names: ['up', 'down'] or ['ud']
    """

    n_sh = len(Sigma_w_sh)

    assert isinstance(spin_names, list)

    with open(filename, 'w') as fo:
        print("# Local self energy at real frequency", file=fo)
        #
        # Column information
        #
        print("# [Column] Data", file=fo)
        print("# [1] Frequency", file=fo)
        icol = 1
        for ish in range(n_sh):
            for sp in spin_names:
                block_dim = Sigma_w_sh[ish][sp].data.shape[1]
                for iorb, jorb in product(list(range(block_dim)), repeat=2):
                    icol += 1
                    print("# [%d] Re(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, sp, iorb, jorb), file=fo)
                    icol += 1
                    print("# [%d] Im(Sigma_{shell=%d, spin=%s, %d, %d})" % (icol, ish, sp, iorb, jorb), file=fo)

        #
        # Write data
        #
        omega = [x.real for x in Sigma_w_sh[0].mesh]
        for iom in range(len(omega)):
            print("%f " % omega[iom], end="", file=fo)
            for ish in range(n_sh):
                for sp in spin_names:
                    block_dim = Sigma_w_sh[ish][sp].data.shape[1]
                    for iorb, jorb in product(list(range(block_dim)), repeat=2):
                        print("%f %f " % (Sigma_w_sh[ish][sp].data[iom, iorb, jorb].real,
                                          Sigma_w_sh[ish][sp].data[iom, iorb, jorb].imag), end="", file=fo)
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


def load_Sigma_iw_sh_txt(filename, Sigma_iw_sh, spin_names, atol_omega=1e-2, interpolation=True):
    """
    Load a list of self-energy from a text file.
    Tails will be set to zero.

    :param filename:
    :param Sigma_iw_sh: All elements must be allocated to correct shapes
    :param spin_names: ['up', 'down'] or ['ud']
    :param atol: absolute tolerance for imaginary frequencies
    :param interpolation: If true and stored data are for a different beta, the data are interpolated along beta axis.
    """

    n_sh = len(Sigma_iw_sh)

    assert isinstance(spin_names, list)

    data = numpy.loadtxt(filename)

    omega_in = data[:, 0]
    omega_out = numpy.array([complex(x) for x in Sigma_iw_sh[0].mesh]).imag
    nomega_out = len(omega_out)

    if not interpolation:
        if not numpy.allclose(omega_in, omega_out, rtol=1, atol=atol_omega):
            raise RuntimeError("Mesh is not compatible!")

    from scipy import interpolate
    interp = interpolate.interp1d(omega_in, data[:,1:], kind='nearest',
        axis=0, fill_value=(data[0,1:], data[-1,1:]), bounds_error=False)

    for iom in range(nomega_out):
        data_interp = interp(omega_out[iom])
        icol = 0
        for ish in range(n_sh):
            for sp in spin_names:
                block_dim = Sigma_iw_sh[ish][sp].data.shape[1]
                for iorb, jorb in product(list(range(block_dim)), repeat=2):
                    re = data_interp[icol]
                    icol += 1
                    imag = data_interp[icol]
                    icol += 1
                    Sigma_iw_sh[ish][sp].data[iom, iorb, jorb] = complex(re, imag)
    if data.shape[1] != \
        numpy.sum([2*Sigma_iw_sh[ish][sp].data.shape[1]**2 for ish in range(n_sh) for sp in spin_names]) + 1:
        raise RuntimeError("Dimensions of the input data are wrong!")

def complex_to_float_array(a):
    a = numpy.ascontiguousarray(a)
    return a.view(float).reshape(a.shape + (2,))

def float_to_complex_array(a):
    a = numpy.ascontiguousarray(a)
    if numpy.iscomplexobj(a):
        return a
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

    #assert type(g) not in [Gf, GfImFreq], 'Type {} is not supported by save_giw'.format(type(g))

    h5file[path + '/__version'] = 'DCore_GfImFreq_v1'
    h5file[path + '/data'] = complex_to_float_array(g.data)
    h5file[path + '/wn'] = numpy.array([complex(x) for x in g.mesh]).imag


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

    assert isinstance(g, Gf), 'Type {} is not supported by save_giw'.format(type(g))
    version_str = h5file[path + '/__version'][()]
    if isinstance(version_str, bytes):
        version_str = version_str.decode('utf-8')
    assert version_str == 'DCore_GfImFreq_v1'

    g.data[...] = float_to_complex_array(h5file[path + '/data'][()])

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
        n_points = g.data.shape[0]//2
        for i in range(n_points):
            diff = numpy.amax(numpy.abs(g.data[i + n_points, :, :]-g.data[n_points - i - 1, :, :].conj().transpose()))
            max_diff = max(max_diff, diff)
            if not check_only:
                g.data[i + n_points, :, :] = g.data[n_points - i - 1, :, :].conj().transpose()
    return max_diff


def _to_numpy_array(g):
    """
    Convert BlockGf object to numpy.
    If there are two blocks, we assume that they are up and down components.
    """

    block_names = gf_block_names(g.n_blocks==1)
    block_sizes = [get_block_size(g[name]) for name in block_names]
    n_spin_orbital = numpy.sum(block_sizes)

    # FIXME: Bit ugly
    n_data = g[block_names[0]].data.shape[0]

    data = numpy.zeros((n_data, n_spin_orbital, n_spin_orbital), dtype=complex)
    offset = 0
    for ib, name in enumerate(block_names):
        block = g[name]
        block_dim = block_sizes[ib]
        data[:, offset:offset + block_dim, offset:offset + block_dim] = block.data
        offset += block_dim
    return data


def _assign_from_numpy_array(g, data):
    """
    Does inversion of to_numpy_array
    """

    if g.n_blocks > 2:
        raise RuntimeError("n_blocks must be 1 or 2.")

    block_names = gf_block_names(g.n_blocks==1)
    n_spin_orbital = numpy.sum([len(block.indices) for name, block in g])

    assert data.shape[0] == g[block_names[0]].data.shape[0]

    norb = n_spin_orbital//2

    offset = 0
    for name in block_names:
        block = g[name]
        block_dim = len(block.indices)
        block.data[:,:,:] = data[:, offset:offset + block_dim, offset:offset + block_dim]
        for i in range(block.data.shape[0]):
            block.data[i, :, :] = 0.5 * (block.data[i, :, :] + block.data[i, :, :].transpose().conj())
        offset += block_dim


def _symmetrize(Sigma_iw, s):
    Sigma_iw_symm = Sigma_iw.copy()
    for bname, gf in Sigma_iw:
        U = s[bname].copy()
        dim_symm = 1
        gf_symm = gf.copy()
        while True:
            gf_rot = gf.copy()
            gf_rot.from_L_G_R(U.transpose().conjugate(), gf, U)
            gf_symm += gf_rot
            dim_symm += 1
            U = numpy.dot(U, s[bname])
            if numpy.allclose(U, numpy.identity(U.shape[0])):
                break
        gf_symm /= dim_symm
        Sigma_iw_symm[bname] = gf_symm
    return Sigma_iw_symm

def symmetrize(Sigma_iw, generators):
    """
    Symmetrize Green's function using generators

    Parameters
    ----------
    Sigma_iw: BlockGf
        self-energy
    generators: list of numpy 2D array
        generators of symmetrization

    Returns
    -------
        Symmetrized self-energy.

    """
    assert isinstance(generators, list)

    Sigma_iw_symm = Sigma_iw.copy()

    for s in generators:
        Sigma_iw_symm = _symmetrize(Sigma_iw_symm, s)

    return Sigma_iw_symm

def mpi_split(work_size, comm_size):
    """
    Make Sigma(iw_n) or G(iwn_n) hermite
    Return max difference.
    """
    base = work_size // comm_size
    leftover = int(work_size % comm_size)

    sizes = numpy.ones(comm_size, dtype=int) * base
    sizes[:leftover] += 1

    offsets = numpy.zeros(comm_size, dtype=int)
    offsets[1:] = numpy.cumsum(sizes)[:-1]

    return sizes, offsets

def expand_path(exec_path):
    """
    Expand relative path and command into full path

    Parameters
    ----------
    exec_path: str
        path or command

    Returns
    -------
    full_path: str
        Full path of exec_path

    """

    full_path = os.path.expandvars(exec_path)  # expand environment variables
    full_path = shutil.which(full_path)  # return full path
    if full_path is None:
        print(f"ERROR: {exec_path} does not exist. Set exec_path properly!", file=sys.stderr)
        sys.exit(1)

    return full_path


def _calc_density(gf):
    """Calculte density by Matsubara summation

    Parameters
    ----------
    gf : GfImFreq

    Returns
    -------
    total_density: float

    """
    assert isinstance(gf, Gf)

    beta = gf.mesh.beta

    # sum over iw and trace over orb of g.data[iw, orb1, orb2]
    density = numpy.sum(numpy.trace(gf.data, axis1=1, axis2=2)) / beta

    # contribution of 1/iw
    density += 0.5 * gf.data.shape[1]

    # equivalent to
    # iw = numpy.array([g.mesh(n) for n in range(g.mesh.first_index(), g.mesh.last_index()+1)])
    # assert iw.shape[0] == g.data.shape[0]
    # density += (0.5 - numpy.sum(1 / iw) / beta) * g.data.shape[1]

    return density


def calc_total_density(block_gf):
    """Calculte total density by Matsubara summation

    Parameters
    ----------
    block_gf : BlockGf

    Returns
    -------
    total_density: float

    """
    assert isinstance(block_gf, BlockGf)

    return numpy.sum([_calc_density(g) for _, g in block_gf])


def _calc_density_matrix(gf):
    """Calculte density matrix by Matsubara summation

    Parameters
    ----------
    gf : GfImFreq

    Returns
    -------
    density_matrix: numpy.ndarray

    """
    assert isinstance(gf, Gf)

    beta = gf.mesh.beta

    # sum over iw of g.data[iw, orb1, orb2]
    dm = numpy.sum(gf.data, axis=0) / beta

    # contribution of 1/iw
    dm += numpy.identity(gf.data.shape[1]) * 0.5

    return dm


def calc_density_matrix(block_gf):
    """Calculte density matrix by Matsubara summation

    Parameters
    ----------
    block_gf : BlockGf

    Returns
    -------
    density_matrix: dict(numpy.ndarray)

    """
    assert isinstance(block_gf, BlockGf)

    return {bname: _calc_density_matrix(g) for bname, g in block_gf}
