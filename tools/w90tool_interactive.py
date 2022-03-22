from __future__ import print_function
import numpy
from itertools import product
import datetime

try:
    input = raw_input  # for python2
except NameError:
    pass


help_commands = """
------------------------------------------------------------------
List of available commands.

(c)ommands
    show command list

(i)nfo
    show info

(h)amiltonian rx ry rz
    show hamiltonian H(rx, ry, rz)

(do)uble spin-major|orbital-major
    double the bases
      spin-major    : 0 1 2 ... -> 0 1 2 ... 0 1 2 ...
      orbital-major : 0 1 2 ... -> 0 0 1 1 2 2 ...

(r)earrange 4 5 2 3
    rearrange bases
      0 1 2 3 4 5 ... -> 4 5 2 3
    Each integer should be >=0 and <NWann, and no duplicate.
    If integers provided are fewer than NWann, some bases will be truncated.

(t)ransform 0-3 file_unitary_matrix
    transform bases
      H_{i,j}(R) -> W^{dag}_{i,i'} H_{i',j'}(R) W_{j',j}
    0-3 indicates bases from 0 to 3 (4 bases) will be transformed.
    The file contains a unitary matrix W (4x4 in this example).

(p)otential 0-3 file_potential
    add a local potential
    The file contains a hermitian matrix H_{pot} (4x4 in this example).

(di)agonalize 0-3:
    transform the bases so that the 0-3 block of H(R=0) becomes diagonal.

(s)ave filename
    save the resultant hamiltonian into a file

(q)uit
------------------------------------------------------------------
"""


def is_unitary(matrix):
    assert matrix.shape[0] == matrix.shape[1]
    # M^dag.M = 1
    return numpy.allclose(numpy.dot(matrix.conjugate().T, matrix), numpy.identity(matrix.shape[0]))


def is_hermitian(matrix):
    assert matrix.shape[0] == matrix.shape[1]
    # M^dag = M
    return numpy.allclose(matrix.conjugate().T, matrix)


class Wannier90( object ):
    def __init__(self, file_name):
        """
        """

        with open(file_name, 'r') as f:
            f.readline()
            self.Nwann = int(f.readline())
            self.nrpts = int(f.readline())

            num_lines = self.nrpts // 15
            if self.nrpts%15 != 0:
                num_lines += 1
            ndgen = []
            for iline in range(num_lines):
                ndgen.extend(f.readline().split())
            self.ndgen = numpy.array(ndgen, dtype=int)

            self.HamR = numpy.zeros((self.nrpts, self.Nwann, self.Nwann), dtype=complex)
            self.irvec = numpy.zeros((self.nrpts, 3), dtype=int)
            for ir in range(self.nrpts):
                for j, i in product(range(self.Nwann), repeat=2):
                    i1, i2, i3, i_in, j_in, hr_real, hr_imag = f.readline().split()
                    if i==0 and j==0:
                        self.irvec[ir,0] = int(i1)
                        self.irvec[ir,1] = int(i2)
                        self.irvec[ir,2] = int(i3)
                    else:
                        assert self.irvec[ir,0] == int(i1)
                        assert self.irvec[ir,1] == int(i2)
                        assert self.irvec[ir,2] == int(i3)
                    assert i == int(i_in)-1
                    assert j == int(j_in)-1
                    self.HamR[ir, i, j] = complex(float(hr_real), float(hr_imag))

    def info(self):
        print("Num of Wannier functions = ", self.Nwann)
        print("Num of R points = ", self.nrpts)

    def print(self, i1, i2, i3):
        rvec = numpy.array([i1, i2, i3], dtype=int)
        ir = self.get_ir(rvec)
        h_r = self.HamR[ir]
        print(h_r)

    def get_ir(self, rvec):
        for ir in range(self.nrpts):
            if (self.irvec[ir] == rvec).all():
                return ir

        raise ValueError("coordinate {} {} {} was not found".format(i1, i2, i3))

    def get_Hk(self, kvec):
        """
        Compute H(k)
        :param kvec: (float, float, float). Fraction coordinates in k space.
        :return: matrix of H(k)
        """

        Hk = numpy.zeros((2*self.norb, 2*self.norb),dtype=complex)
        for iR in range(self.nrpts):
            factor = numpy.exp(2J*numpy.pi*(self.irvec[iR,0]*kvec[0]+self.irvec[iR,1]*kvec[1]+self.irvec[iR,2]*kvec[2]))
            Hk += self.HamR_full[iR,:,:] * factor / self.ndgen[iR]
        return Hk

    def save(self, filename):
        # f = file(filename, 'w')
        with open(filename, 'w') as f:
            # print('DUMMY', file=f)
            print(datetime.datetime.now(), file=f)
            # print(2*self.norb, file=f)
            print(self.Nwann, file=f)
            print(self.nrpts, file=f)
            for k in range(self.nrpts):
                print(self.ndgen[k], file=f, end=' ')
                if k % 15 == 14 or k == self.nrpts - 1:
                    print('', file=f)

            for k in range(self.nrpts):
                # for j, i in product(range(2 * self.norb), repeat=2):
                for j, i in product(range(self.Nwann), repeat=2):
                    print("{} {} {}  {} {}  {} {}".format(
                        self.irvec[k, 0], self.irvec[k, 1], self.irvec[k, 2],
                        i + 1, j + 1,
                        self.HamR[k, i, j].real, self.HamR[k, i, j].imag), file=f)

    def _check_index(self, index):
        if index >= self.Nwann:
            raise ValueError("Index {} is out of range".format(index))

    def _check_indices(self, indices):
        for index in indices:
            self._check_index(index)

    def double(self, spin_orbital):
        if spin_orbital == 'spin-major':
            ham_new = numpy.zeros((self.nrpts, 2, self.Nwann, 2, self.Nwann), dtype=complex)
            for isp in range(2):
                ham_new[:, isp, :, isp, :] = self.HamR
        elif spin_orbital == 'orbital-major':
            ham_new = numpy.zeros((self.nrpts, self.Nwann, 2, self.Nwann, 2), dtype=complex)
            for isp in range(2):
                ham_new[:, :, isp, :, isp] = self.HamR

        self.Nwann *= 2
        self.HamR = ham_new.reshape((self.nrpts, self.Nwann, self.Nwann))

    def rearrange(self, indices):
        self._check_indices(indices)

        # check duplicate
        if len(indices) != len(set(indices)):
            raise ValueError("Some indices duplicate")

        # if max(indices) >= self.Nwann:
        #     raise ValueError("Some indices out of range")

        nwann = len(indices)

        if nwann == self.Nwann:
            ham_new = self.HamR[:, indices, :][:, :, indices]
        else:
            ham_new = numpy.zeros((self.nrpts, nwann, nwann), dtype=complex)
            for ir in range(self.nrpts):
                for j, i in product(range(nwann), repeat=2):
                    ham_new[ir, j, i] = self.HamR[ir, indices[j], indices[i]]

        self.Nwann = nwann
        self.HamR = ham_new

    def transform(self, w_matrix, istart, iend):
        """
        H_{i,j}(R) -> W^{dag}_{i,i'} H_{i',j'}(R) W_{j',j}

        W: unitary matrix
        """
        self._check_indices([istart, iend])

        # check if 'unitary_matrix' is unitary
        if not is_unitary(w_matrix):
            raise ValueError("The transformation matrix is not unitary")

        mat_full = numpy.identity(self.Nwann, dtype=complex)
        mat_full[istart:iend+1, istart:iend+1] = w_matrix

        for ir in range(self.nrpts):
            h = self.HamR[ir]
            h_transformed = numpy.dot(mat_full.conjugate().T, numpy.dot(h, mat_full))
            self.HamR[ir] = h_transformed

    def add_potential(self, v_matrix, istart, iend):
        """
        H_{i,j}(R=) += V_{i,j}

        V: hermitian matrix
        """
        self._check_indices([istart, iend])

        # check if 'potential_matrix' is hermitian
        if not is_hermitian(v_matrix):
            raise ValueError("The potential matrix is not hermitian")

        ir = self.get_ir(numpy.array([0, 0, 0]))
        self.HamR[ir, istart:iend+1, istart:iend+1] += v_matrix

    def diagonalize(self, istart, iend):
        self._check_indices([istart, iend])
        ir = self.get_ir(numpy.array([0, 0, 0]))
        h0 = self.HamR[ir, istart:iend+1, istart:iend+1]
        # print(h0.shape)

        # TODO: check eigvec U^dag.H.U or U.H.U^dag

        eigval, eigvec = numpy.linalg.eigh(h0)
        print("Eigenvalues =", eigval)

        self.transform(eigvec, istart, iend)


def _read_local_coordinate(filename):
    f = open(filename, 'r')

    def _readline():
        while True:
            line = f.readline()
            if not line or not line.startswith('#'):
                break
        return line

    ncor_sh = int(_readline())
    for ish in range(ncor_sh):
        print("ish = {}".format(ish))
        x_axis = numpy.array(map(float, _readline().split()))
        y_axis = numpy.array(map(float, _readline().split()))
        z_axis = numpy.array(map(float, _readline().split()))
        print(" x_axis = ", x_axis)
        print(" y_axis = ", y_axis)
        print(" z_axis = ", z_axis)

        x_axis /= numpy.linalg.norm(x_axis)
        y_axis /= numpy.linalg.norm(y_axis)
        z_axis /= numpy.linalg.norm(z_axis)

        u_mat = numpy.array([x_axis, y_axis, z_axis]).transpose()

        if not numpy.allclose(u_mat, u_mat.transpose()):
            raise RuntimeError("Local coordinate axes are not orthogonal!")
    f.close()


def _check_num_args(commands, num):
    if len(commands) - 1 < num:
        raise ValueError("Command '{}' requires {} argument(s)".format(commands[0], num))


def _get_istart_iend(arg):
    """Get two integers from hyphen-separated numbers, e.g., 0-3 --> 0, 3"""
    range_split = arg.split('-')
    assert len(range_split) == 2, "Failed to split '{}' into to integers".format(arg)
    istart = int(range_split[0])
    iend = int(range_split[1])
    assert istart < iend, "{} < {}".format(istart, iend)
    return istart, iend


def _read_file(filename, shape):
    print("Reading file '{}'".format(filename))
    matrix = numpy.loadtxt(filename).view(complex)
    # print(matrix.shape)
    assert matrix.shape == shape, \
        "The matrix size {} in file '{}' does not match {}".format(matrix.shape, filename, shape)
    return matrix


if __name__ == '__main__':
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(
        prog='w90tool.py',
        description='An interactive tool for manipulating a wannier90 format file' + help_commands,
        # usage='$ python w90tool.py seedname',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument('filename', action='store', default=None, type=str, help="Input filename. Normally, *seedname*_hr.dat")
    parser.add_argument('--terminate-if-exception', '-t', action='store_true', help="Terminate the program when an exception occurs. This flag is recommended to activate when commands are given from a file.")

    args = parser.parse_args()
    # print(args)

    filename = args.filename
    if os.path.isfile(filename) is False:
        print("file '{}' not found".format(filename))
        sys.exit(-1)

    w90 = Wannier90(args.filename)
    w90.info()

    print(help_commands, end="")

    while True:
        try:
            print("\ncommand? ", end="")
            commands_str = input()
            # print(commands_str)
            commands = commands_str.split(' ')
            # print(commands)

            if commands[0] in ('info', 'i'):
                _check_num_args(commands, 0)
                w90.info()

            elif commands[0] in ('comamnds', 'c'):
                _check_num_args(commands, 0)
                print(help_commands, end="")

            elif commands[0] in ('hamiltonian', 'h'):
                _check_num_args(commands, 3)
                i1 = int(commands[1])
                i2 = int(commands[2])
                i3 = int(commands[3])
                w90.print(i1, i2, i3)

            elif commands[0] in ('double', 'do'):
                _check_num_args(commands, 1)
                spin_orbital = commands[1]
                if spin_orbital not in ('spin-major', 'orbital-major'):
                    raise ValueError("'{}' should be either 'spin-major' or 'orbital-major'".format(spin_orbital))
                w90.double(spin_orbital)
                print("Succeeded")

            elif commands[0] in ('rearrange', 'r'):
                _check_num_args(commands, 1)

                indices = numpy.array(commands[1:], dtype=int)
                print("indices = {}".format(indices))

                w90.rearrange(indices)
                print("Succeeded")

            elif commands[0] in ('transform', 't'):
                _check_num_args(commands, 2)

                istart, iend = _get_istart_iend(commands[1])
                print("Bases from {} to {} will be transformed".format(istart, iend))

                file_basis = commands[2]
                matrix = _read_file(file_basis, (iend - istart + 1, iend - istart + 1))
                # TODO: format of the matrix U, U^T, or U^dag?

                w90.transform(matrix, istart, iend)
                print("Succeeded")

            elif commands[0] in ('potential', 'p'):
                _check_num_args(commands, 2)

                istart, iend = _get_istart_iend(commands[1])
                print("Add a potential to H(R=0) bases from {} to {}".format(istart, iend))

                file_pot = commands[2]
                matrix = _read_file(file_pot, (iend - istart + 1, iend - istart + 1))

                w90.add_potential(matrix, istart, iend)
                print("Succeeded")

            elif commands[0] in ('diagonalize', 'di'):
                _check_num_args(commands, 1)

                istart, iend = _get_istart_iend(commands[1])
                print("Diagonalize H(R=0) bases from {} to {}".format(istart, iend))

                w90.diagonalize(istart, iend)
                print("Succeeded")

            elif commands[0] in ('save', 's'):
                _check_num_args(commands, 1)

                fileout = commands[1]
                print("Saving to file '{}'".format(fileout))
                w90.save(fileout)
                print("Succeeded")

            elif commands[0] in ('quit', 'q'):
                _check_num_args(commands, 0)
                break

            else:
                raise ValueError("Command '{}' not supported".format(commands[0]))

        except ValueError as e:
            print("ERROR:", e)
            if args.terminate_if_exception:
                sys.exit(1)

        except AssertionError as e:
            print("AssertionError:", e)
            if args.terminate_if_exception:
                sys.exit(1)

        except IOError as e:
            print("IOError:", e)
            if args.terminate_if_exception:
                sys.exit(1)

        except EOFError:
            print("EOF")
            break
