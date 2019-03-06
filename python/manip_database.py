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
    data_new['n_orbitals'] = numpy.zeros((nk, 1), dtype=int)
    data_new['n_orbitals'][:, 0] = data['n_orbitals'][:, sp[0]] + data['n_orbitals'][:, sp[1]]

    max_n_orbitals = data_new['n_orbitals'].max()
    print("max_n_orbitals =", max_n_orbitals)

    data_new['proj_mat'] = numpy.zeros((nk, 1, n_corr_shells, max_n_orbitals, max_n_orbitals), dtype=complex)
    for k in range(nk):
        offset = 0
        for s in sp:
            norb = data['n_orbitals'][k, s]
            data_new['proj_mat'][k, 0, :, offset:offset+norb, offset:offset+norb]\
                = data['proj_mat'][k, s, :, 0:norb, 0:norb]
            offset += norb

    data_new['hopping'] = numpy.zeros((nk, 1, max_n_orbitals, max_n_orbitals), dtype=complex)
    for k in range(nk):
        offset = 0
        for s in sp:
            norb = data['n_orbitals'][k, s]
            data_new['hopping'][k, 0, offset:offset+norb, offset:offset+norb]\
                = data['hopping'][k, s, 0:norb, 0:norb]
            offset += norb

    for key in keys:
        print("{:11s}:".format(key), data_new[key].shape)

    #
    # set flags
    #
    data_new['SO'] = 1
    data_new['SP'] = 1

    # corr_shells
    corr_shells = data['corr_shells']
    print(corr_shells)
    for crsh in corr_shells:
        assert crsh['SO'] == 0
        crsh['SO'] = 1
    data_new['corr_shells'] = corr_shells

    #
    # save data
    #
    if h5_file_in == h5_file_out:
        data_save = data_new
    else:
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
