

import numpy
from itertools import product

import argparse


def run():
    # python mk_init_se.py --x=1.0 --y=1.0 --z=1.0 --output=init_se.txt --norb=3
    # This will generate an initial self-energy for pointing local spins to (1, 1, 1).
    # Only spin_orbit case is supported.
    parser = argparse.ArgumentParser(
        prog='mk_init_se.py',
        description='make inital guess for self-energy',
        usage='$python mk_init_se.py --x=1.0 --y=1.0 --z=1.0 --output=init_se.txt --norb=3',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('--x', default=0.0, help='x', required=True)
    parser.add_argument('--y', default=0.0, help='y', required=True)
    parser.add_argument('--z', default=0.0, help='z', required=True)
    parser.add_argument('--norb', default=1, help='number of orbitals', required=True)
    parser.add_argument('--output', default='', help='output file', required=True)

    args = parser.parse_args()
    
    Sx = numpy.array([[0, 1], [1,0]], dtype=complex)
    Sy = numpy.array([[0, -1J], [1J, 0]], dtype=complex)
    Sz = numpy.array([[1, 0], [0, -1]], dtype=complex)
    
    x, y, z = float(args.x), float(args.y), float(args.z)
    norb = int(args.norb)

    S = x*Sx + y*Sy + z*Sz

    with open(args.output, 'w') as f:
        for i, (isp, iorb) in enumerate(product(list(range(2)), list(range(norb)))):
            for j, (jsp, jorb) in enumerate(product(list(range(2)), list(range(norb)))):
                if iorb==jorb:
                    print(0, i, j, -S[isp, jsp].real, -S[isp, jsp].imag, file=f)
                else:
                    print(0, i, j, 0.0, 0.0, file=f)
