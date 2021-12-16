
import numpy
from itertools import product

from numpy.testing import assert_allclose
#import sympy
#import pytransform3d
#import re

class Wannier90( object ):
    def __init__(self, file_name, spin_orbital_order='up_up_down_down', verbose=0):
        """
        spin_orbital_order = 'up_down_up_down', 'up_up_down_down' or 'up_up_up_up'.
        For spin_orbital_order == (up_up_up_up), only the up spin elements are given in the input file
        and the spin-full Hamiltonian will be constructed.
        Internally, the ordering of spins and orbitals in HamR is (orb1, up), (orb2, up) ... (orb1, down) ...
        """
        f = open(file_name)
        f.readline()
        self.Nwann = int(f.readline())
        self.nrpts = int(f.readline())

        if verbose > 0:
            print("Num of Wannier functions = ", self.Nwann)
            print("Num of R points = ", self.nrpts)

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
            for j in range(self.Nwann):
                for i in range(self.Nwann):
                    i1, i2, i3, i_in, j_in, hr_real,hr_imag = f.readline().split()
                    if i==0 and j==0:
                        self.irvec[ir,0] = i1
                        self.irvec[ir,1] = i2
                        self.irvec[ir,2] = i3
                    assert i == int(i_in)-1
                    assert j == int(j_in)-1
                    self.HamR[ir, i, j] = complex(float(hr_real), float(hr_imag))

        if spin_orbital_order == 'up_down_up_down':
            # (R, orb, spin, orb, spin) => (R, spin, orb, spin, orb)
            self.HamR_full = self.HamR.reshape((self.nrpts, self.Nwann//2, 2, self.Nwann//2, 2)).transpose((0, 2, 1, 4, 3))
            self.HamR_full = self.HamR_full.reshape((self.nrpts, self.Nwann, self.Nwann))
            self.norb = self.Nwann//2
        elif spin_orbital_order == 'up_up_down_down':
            self.HamR_full = self.HamR
            self.norb = self.Nwann//2
        elif spin_orbital_order == 'up_up_up_up':
            self.norb = self.Nwann
            self.HamR_full = numpy.zeros((self.nrpts, 2, self.norb, 2, self.norb), dtype=complex)
            for isp in range(2):
                self.HamR_full[:, isp, :, isp, :] = self.HamR
            self.HamR_full = self.HamR_full.reshape((self.nrpts, 2 * self.norb, 2 * self.norb))
        else:
            raise RuntimeError("Invalid spin_orbital_order: {}".format(spin_orbital_order))

    def get_Hk(self, kvec, check_hermite=False):
        """
        Compute H(k)
        :param kvec: (float, float, float). Fraction coordinates in k space.
        :return: matrix of H(k)
        """

        Hk = numpy.zeros((2*self.norb, 2*self.norb),dtype=complex)
        for iR in range(self.nrpts):
            factor = numpy.exp(2J*numpy.pi*(self.irvec[iR,0]*kvec[0]+self.irvec[iR,1]*kvec[1]+self.irvec[iR,2]*kvec[2]))
            Hk += self.HamR_full[iR,:,:] * factor / self.ndgen[iR]
        if check_hermite:
            assert_allclose(Hk, Hk.T.conjugate())
        return Hk

    def save(self, filename):
        f = open(filename, 'w')
        print('DUMMY', file=f)
        print(2*self.norb, file=f)
        print(self.nrpts, file=f)
        for k in range(self.nrpts):
            print(self.ndgen[k], file=f, end=' ')
            if k % 15 == 14 or k == self.nrpts - 1:
                print('', file=f)

        for k in range(self.nrpts):
            for j, i in product(list(range(2 * self.norb)), repeat=2):
                print("{} {} {}  {} {}  {} {}".format(
                    self.irvec[k, 0], self.irvec[k, 1], self.irvec[k, 2],
                    i + 1, j + 1,
                    self.HamR_full[k, i, j].real, self.HamR_full[k, i, j].imag), file=f)

        f.close()

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
        x_axis = numpy.array(list(map(float, _readline().split())))
        y_axis = numpy.array(list(map(float, _readline().split())))
        z_axis = numpy.array(list(map(float, _readline().split())))
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


def run():
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(
        prog='w90tool.py',
        description='.',
        usage='$ python w90tool.py seedname',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument('seedname',
                        action='store',
                        default=None,
                        type=str,
                        help="seedname")
    parser.add_argument('--inputorder', default=None, help='up_down_up_down, up_up_down_down or up_up_up_up', required=True)
    #parser.add_argument('--localaxis', default=None, help='Name of file containing local coordinate axes', required=True)
    parser.add_argument('--output', default=None, help='Output file name', required=False)

    args = parser.parse_args()

    if os.path.isfile(args.seedname+'_hr.dat') is False:
        print("{}_hr.dat is not found".format(args.seedname))
        sys.exit(-1)

    w90 = Wannier90(args.seedname+'_hr.dat', spin_orbital_order=args.inputorder)

    #_read_local_coordinate(args.localaxis)

    if not args.output is None:
        w90.save(args.output)
