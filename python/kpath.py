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

import numpy
from collections import namedtuple

from .program_options import parse_knode, parse_bvec


# x-position of k-node (symmetric point of k-vector)
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
