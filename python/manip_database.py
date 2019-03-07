#!/usr/bin/env python
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
import os
from pytriqs.archive.hdf_archive import HDFArchive


def expand_block_diag_matrix(mat1, mat2):
    """
    make a large-size matrix combining 2 matrices as diag(mat1, mat2)
    """
    dim1_l = mat1.shape[0]
    dim1_r = mat1.shape[1]
    dim2_l = mat2.shape[0]
    dim2_r = mat2.shape[1]

    mat = numpy.zeros((dim1_l + dim2_l, dim1_r + dim2_r), dtype=complex)
    mat[0: dim1_l, 0: dim1_r] = mat1[:, :]
    mat[dim1_l: dim1_l + dim2_l, dim1_r: dim1_r + dim2_r] = mat2[:, :]
    return mat


def double_matrix(mat):
    """
    make a double-size block matrix, diag(mat, mat)
    """
    return expand_block_diag_matrix(mat, mat)


def turn_on_spin_orbit(h5_file_in, h5_file_out):
    if not os.path.exists(h5_file_in):
        print("file '{}' not found".format(h5_file_in))
        exit(1)

    #
    # read data
    #
    with HDFArchive(h5_file_in, 'r') as ar:
        ar = ar['dft_input']
        data = {key: ar[key] for key in ar.keys()}
    print(data.keys())

    n_corr_shells = data['n_corr_shells']
    nk = data['n_k']

    # check flags
    SP = data['SP']  # spin polarized
    SO = data['SO']  # spin-orbit
    print("(SP, SO) = ({}, {})".format(SP, SO))
    assert (SP, SO) == (0, 0) or (SP, SO) == (1, 0)

    sp = [0, 0] if (SP, SO) == (0, 0) else [0, 1]

    #
    # expand arrays
    #
    keys = ('n_orbitals', 'proj_mat', 'hopping')  # arrays to be updated
    for key in keys:
        print("{:11s}:".format(key), data[key].shape)

    data_new = {}

    # n_orbitals
    data_new['n_orbitals'] = numpy.zeros((nk, 1), dtype=int)
    data_new['n_orbitals'][:, 0] = data['n_orbitals'][:, sp[0]] + data['n_orbitals'][:, sp[1]]

    max_n_orbitals = data_new['n_orbitals'].max()
    print("max_n_orbitals =", max_n_orbitals)

    max_corr_shell_dim = max([ corr_shell['dim'] for corr_shell in data['corr_shells'] ])
    max_corr_shell_dim *= 2
    print("max_corr_shell_dim =", max_corr_shell_dim)

    # proj_mat
    data_new['proj_mat'] = numpy.zeros((nk, 1, n_corr_shells, max_corr_shell_dim, max_n_orbitals), dtype=complex)
    for k in range(nk):
        for j in range(n_corr_shells):
            data_new['proj_mat'][k, 0, j] \
                = expand_block_diag_matrix(data['proj_mat'][k, sp[0], j], data['proj_mat'][k, sp[1], j])

    # hopping
    data_new['hopping'] = numpy.zeros((nk, 1, max_n_orbitals, max_n_orbitals), dtype=complex)
    for k in range(nk):
        data_new['hopping'][k, 0] \
            = expand_block_diag_matrix(data['hopping'][k, sp[0]], data['hopping'][k, sp[1]])

    for key in keys:
        print("{:11s}:".format(key), data_new[key].shape)

    # rot_mat
    data_new['rot_mat'] = [double_matrix(rot) for rot in data['rot_mat']]

    # T
    data_new['T'] = [double_matrix(t) for t in data['T']]

    #
    # set flags
    #
    data_new['SO'] = 1
    data_new['SP'] = 1

    # corr_shells
    corr_shells = data['corr_shells']
    # print(corr_shells)
    for crsh in corr_shells:
        assert crsh['SO'] == 0
        crsh['SO'] = 1
        crsh['dim'] *= 2
    data_new['corr_shells'] = corr_shells
    # print(corr_shells)

    #
    # print changes in database
    #
    def print_formated(obj):
        if isinstance(obj, numpy.ndarray):
            print("    ", obj.shape)
        elif isinstance(obj, list):
            for x in obj:
                print_formated(x)
        else:
            print("    ", obj)

    print("Summary of changes in the database")
    for key, d in data_new.items():
        print(" -", key)
        print_formated(data[key])
        print_formated(d)

    #
    # save data
    #
    if h5_file_in == h5_file_out:
        # overwrite only updated components
        data_save = data_new
    else:
        # write all components
        data.update(data_new)
        data_save = data

    with HDFArchive(h5_file_out, 'a') as ar:
        if 'dft_input' not in ar:
            ar.create_group('dft_input')
        ar = ar['dft_input']

        for key, d in data_save.items():
            ar[key] = d


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_h5_file', help='input hdf5 file')
    parser.add_argument('-o', help='output hdf5 file. If not given, input data are overwritten')
    args = parser.parse_args()
    print(args)

    if args.o is None:
        output_h5_file = args.input_h5_file
    else:
        output_h5_file = args.o

    turn_on_spin_orbit(args.input_h5_file, output_h5_file)
