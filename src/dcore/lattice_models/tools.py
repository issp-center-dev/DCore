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


import numpy
import sys
from collections import namedtuple

from dcore._dispatcher import HDFArchive

from ..tools import pauli_matrix
from ..program_options import parse_knode, parse_bvec


def _drop_small_vals(z, eps=1e-10):
    z_real = z.real if numpy.abs(z.real) > eps else 0.0
    z_imag = z.imag if numpy.abs(z.imag) > eps else 0.0
    return complex(z_real, z_imag)


def print_spin_orbital_matrix(mat, file, print_offset=0):
    norb = mat.shape[1]

    for isp in range(2):
        for iorb in range(norb):
            print(' '*print_offset, end='')
            for jsp in range(2):
                for jorb in range(norb):
                    z = _drop_small_vals(mat[isp, iorb, jsp, jorb])
                    print('({0:>9.2e},{1:>9.2e})'.format(z.real, z.imag), end=' ', file=file)
                if jsp==0:
                    print('  ', file=file, end='')
            print('', file=file)
        print('', file=file) 


def print_local_fields(h5_file, corr_shell_dims=None, subgrp='dft_input'):
    """
    Print local fields of H(R=0)
    :param h5_file: input file for DFTTools
    :param subgrp:
    """

    with HDFArchive(h5_file, 'r') as f:
        Hk = f[subgrp]['hopping'][()]
        SO = f[subgrp]['SO']
        SP = f[subgrp]['SP']
        bz_weights = f[subgrp]['bz_weights'][()]
        n_corr_sh = f[subgrp]['n_corr_shells']
        dims_corr_sh = numpy.array([f[subgrp]['corr_shells'][ish]['dim'] for ish in range(n_corr_sh)])

    if (SO==1 and SP==0) or (SO==0 and SP==1):
        raise RuntimeError("SO={} and SP={} are not supported by DCore!".format(SO, SP))

    nk = Hk.shape[0]
    spin_block_dim = Hk.shape[2]
    if SO==0:
        Hk_ud = numpy.zeros((nk, 2, spin_block_dim, 2, spin_block_dim), dtype=complex)
        for isp in range(2):
            Hk_ud[:, isp, :, isp, :] = Hk[:, 0, :, :]
        Hk_ud = Hk_ud.reshape((nk, 2*spin_block_dim, 2*spin_block_dim))
        norb = spin_block_dim
        num_spin_orb_corr_sh = 2*dims_corr_sh
    else:
        Hk_ud = Hk.reshape((nk, spin_block_dim, spin_block_dim))
        norb = spin_block_dim//2
        num_spin_orb_corr_sh = dims_corr_sh

    # Hk_ud (nk, 2*num_orb, 2*num_orb)
    # Check only H(R=0)
    H0_ud = numpy.einsum('kij,k->ij', Hk_ud, bz_weights)
    H0_ud = H0_ud.reshape((2, norb, 2, norb))

    print('')
    print('---local fields (w/o local potential)')
    pauli_mat = pauli_matrix()
    for iorb in range(norb):
        # alpha = x, y, z
        h = numpy.empty(3)
        for alpha in range(3):
            h[alpha] = 0.5 * numpy.trace(numpy.dot(H0_ud[:, iorb, :, iorb], pauli_mat[alpha])).real
        print("    orbital {} : hx, hy, hz = {} {} {}".format(iorb, *h))

    print('')
    print('---intra shell Hamiltonian (w/o local potential)')
    # (spin, orb, spin, orb) => (orb, spin, orb, spin)
    H0_ud2 = H0_ud.transpose((1, 0, 3, 2)).reshape((2*norb, 2*norb))
    offset = 0
    print('    ordering: up ... up, down ... down')
    for ish in range(n_corr_sh):
        print(' '*4, 'corr_shell=', ish)
        block_size = num_spin_orb_corr_sh[ish]
        H0_block = H0_ud2[offset:offset+block_size, offset:offset+block_size]
        # (orb, spin, orb, spin) => (spin, orb, spin, orb)
        H0_block = H0_block.reshape((block_size//2, 2, block_size//2, 2)).transpose((1, 0, 3, 2))
        print_spin_orbital_matrix(H0_block, sys.stdout, 6)
        evals, evecs = numpy.linalg.eigh(H0_block.reshape((block_size, block_size)))
        print('      eigenvalues: ', evals)
        print('')
        offset += block_size


# k-node information for plot
#    x : float
#    label : string
XNode = namedtuple('XNode', ('x', 'label'))


def gen_kpath(params):
    """
    Generate k-path

    Parameters
    ----------
    params : dict

    Returns
    -------
    xk : numpy.ndarray
        x-mesh for the band plot

    xnode : list of XNode

    kvec : numpy.ndarray
        k-vectors on the k-path

    """

    knode = parse_knode(params["tool"]["knode"])
    bvec = parse_bvec(params["model"]["bvec"])
    nk_line = params["tool"]["nk_line"]
    #
    # Generate k-vectors on the k-path
    #   knode --> kvec
    #
    nnode = len(knode)
    n_k = (nnode - 1) * nk_line + 1
    # print("\n   Total number of k =", str(n_k))
    kvec = numpy.zeros((n_k, 3), numpy.float_)
    ikk = 0
    for inode in range(nnode - 1):
        for ik in range(nk_line + 1):
            if inode != 0 and ik == 0:
                continue
            kvec[ikk] = float((nk_line - ik)) * knode[inode].kvec + float(ik) * knode[inode + 1].kvec
            kvec[ikk] *= 2.0 * numpy.pi / float(nk_line)
            ikk += 1
    #
    # Generate x values on the k-path (horizontal axis in a band plot)
    #   knode --> xk, xnode
    #
    xk = numpy.zeros(n_k, numpy.float_)
    xnode = []
    xk[0] = 0.0
    ikk = 0
    for inode in range(nnode - 1):
        dk = knode[inode + 1].kvec - knode[inode].kvec
        dk_cart = numpy.dot(dk, bvec)
        klength = numpy.sqrt(numpy.dot(dk_cart, dk_cart)) / nk_line
        xnode.append(XNode(x = xk[ikk], label = knode[inode].label))
        for ik in range(nk_line):
            xk[ikk + 1] = xk[ikk] + klength
            ikk += 1
    xnode.append(XNode(x = xk[-1], label = knode[-1].label))

    return xk, xnode, kvec
