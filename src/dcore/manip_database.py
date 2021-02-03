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


import numpy
import os
import copy
from h5 import HDFArchive


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


class H5SpinOrbitOn:
    def __init__(self, h5_file):
        if not os.path.exists(h5_file):
            print("file '{}' not found".format(h5_file))
            exit(1)
        self.h5_file_in = h5_file

        #
        # read data
        #
        self.set_group('dft_input')

        #
        # define values common to all groups
        #
        self.n_corr_shells = self.data['n_corr_shells']
        self.max_corr_shell_dim = max([corr_shell['dim'] for corr_shell in self.data['corr_shells']]) * 2
        print("max_corr_shell_dim =", self.max_corr_shell_dim)

        # store corr_shells for update of shells
        self.store_corr_shells = copy.deepcopy(self.data['corr_shells'])

        # check flags
        SP = self.data['SP']  # spin polarized
        SO = self.data['SO']  # spin-orbit
        print("(SP, SO) = ({}, {})".format(SP, SO))
        assert SO == 0

        self.sp = [0, 0] if SP == 0 else [0, 1]

    def set_group(self, grp):
        """
        return True if succeeded
        """
        if hasattr(self, 'data'):
            del self.data

        with HDFArchive(self.h5_file_in, 'r') as ar:
            if grp not in ar:
                return False
            ar = ar[grp]
            self.data = {key: ar[key] for key in list(ar.keys())}
        print("\n*** grp = " + grp)
        print(list(self.data.keys()))

        if 'n_k' in self.data:
            self.nk = self.data['n_k']
            print("n_k =", self.nk)

        if hasattr(self, '_max_n_orbitals'):
            del self._max_n_orbitals

        self.data_new = {}

        return True

    def save(self, h5_file_out, grp):
        if self.h5_file_in == h5_file_out:
            # overwrite only updated components
            data_save = self.data_new
        else:
            # write all components
            data_save = dict(self.data)
            data_save.update(self.data_new)

        with HDFArchive(h5_file_out, 'a') as ar:
            if grp not in ar:
                ar.create_group(grp)
            ar = ar[grp]

            for key, d in list(data_save.items()):
                ar[key] = d

    def update(self, key):
        assert key in dir(self)

        self.data_new[key] = eval('self.'+key)()

        def print_formated(obj):
            if isinstance(obj, numpy.ndarray):
                print("    ", obj.shape)
            elif isinstance(obj, list):
                for x in obj:
                    print_formated(x)
            else:
                print("    ", obj)

        print(" -", key)
        print_formated(self.data[key])
        print_formated(self.data_new[key])

    @property
    def max_n_orbitals(self):
        if not hasattr(self, '_max_n_orbitals'):
            print("n_orbitals must be set first")
            exit(1)
        return self._max_n_orbitals

    #
    # implementation of updating data
    #
    def n_orbitals(self):
        n_orbitals = numpy.zeros((self.nk, 1), dtype=int)
        n_orbitals[:, 0] = self.data['n_orbitals'][:, self.sp[0]] + self.data['n_orbitals'][:, self.sp[1]]

        self._max_n_orbitals = n_orbitals.max()
        # print("max_n_orbitals =", self._max_n_orbitals)
        return n_orbitals

    def proj_mat(self):
        max_n_orbitals = self.max_n_orbitals
        max_corr_shell_dim = self.max_corr_shell_dim

        proj_mat = numpy.zeros((self.nk, 1, self.n_corr_shells, max_corr_shell_dim, max_n_orbitals), dtype=complex)
        for k in range(self.nk):
            for j in range(self.n_corr_shells):
                proj_mat[k, 0, j] \
                    = expand_block_diag_matrix(self.data['proj_mat'][k, self.sp[0], j], self.data['proj_mat'][k, self.sp[1], j])
        return proj_mat

    def hopping(self):
        max_n_orbitals = self.max_n_orbitals

        hopping = numpy.zeros((self.nk, 1, max_n_orbitals, max_n_orbitals), dtype=complex)
        for k in range(self.nk):
            hopping[k, 0] \
                = expand_block_diag_matrix(self.data['hopping'][k, self.sp[0]], self.data['hopping'][k, self.sp[1]])
        return hopping

    def corr_shells(self):
        corr_shells = copy.deepcopy(self.data['corr_shells'])
        for crsh in corr_shells:
            assert crsh['SO'] == 0
            crsh['SO'] = 1
            crsh['dim'] *= 2
        return corr_shells

    def shells(self):
        atom_l_sort = [(crsh['atom'], crsh['l'], crsh['sort']) for crsh in self.store_corr_shells]

        shells = copy.deepcopy(self.data['shells'])
        for sh in shells:
            if (sh['atom'], sh['l'], sh['sort']) in atom_l_sort:
                sh['dim'] *= 2
        return shells


    def rot_mat(self):
        return [double_matrix(rot) for rot in self.data['rot_mat']]

    def T(self):
        return [double_matrix(t) for t in self.data['T']]

    def SP(self):
        return 1

    def SO(self):
        return 1


def turn_on_spin_orbit(h5_file_in, h5_file_out, update_dft_input=True, update_dft_bands_input=True):

    # 'dft_input'
    h5so = H5SpinOrbitOn(h5_file_in)
    if update_dft_input:
        for key in ['n_orbitals', 'proj_mat', 'hopping', 'corr_shells', 'shells', 'rot_mat', 'T', 'SP', 'SO']:
            h5so.update(key)
        h5so.save(h5_file_out, 'dft_input')

    # 'dft_bands_input'
    if update_dft_bands_input and h5so.set_group('dft_bands_input'):
        for key in ['n_orbitals', 'proj_mat']:
            h5so.update(key)
        h5so.save(h5_file_out, 'dft_bands_input')


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
