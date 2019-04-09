from __future__ import print_function

import numpy
from itertools import product

import argparse

# python mk_init_se.py --mx=1.0 --my=1.0 --mz=1.0 --output=init_se.txt --norb=3
# This will generate an initial self-energy for pointing local spins to (1, 1, 1).
# Only spin_orbit case is supported.
parser = argparse.ArgumentParser(
    prog='mk_init_se.py',
    description='make inital guess for self-energy',
    usage='$python mk_init_se.py --mx=1.0 --my=1.0 --mz=1.0 --output=init_se.txt --norb=3',
    add_help=True,
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument('--mx', default=0.0, help='mx', required=True)
parser.add_argument('--my', default=0.0, help='mx', required=True)
parser.add_argument('--mz', default=0.0, help='mx', required=True)
parser.add_argument('--norb', default=1, help='number of orbitals', required=True)
parser.add_argument('--output', default='', help='output file', required=True)

args = parser.parse_args()

Sx = numpy.array([[0, 1], [1,0]], dtype=complex)
Sy = numpy.array([[0, -1J], [1J, 0]], dtype=complex)
Sz = numpy.array([[1, 0], [0, -1]], dtype=complex)

mx, my, mz = float(args.mx), float(args.my), float(args.mz)
norb = int(args.norb)

S = mx*Sx + my*Sy + mz*Sz

with open(args.output, 'w') as f:
    for i, (isp, iorb) in enumerate(product(range(2), range(norb))):
        for j, (jsp, jorb) in enumerate(product(range(2), range(norb))):
            if iorb==jorb:
                print(0, i, j, -S[isp, jsp].real, -S[isp, jsp].imag, file=f)
            else:
                print(0, i, j, 0.0, 0.0, file=f)
