#!/usr/bin/python
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
import sys
import numpy


def res2wan(name_in, name_out):
    #
    print("  Convert from \"{0}\" tp \"{1}\"".format(name_in, name_out))
    #
    # Input
    #
    with open(name_in, 'r') as f:
        for ii in range(3):
            line = f.readline()  # Skip
            print("    "+line, end="")
        line = f.readlines()  # Skip
    #
    # Parse from RESPACK format
    #
    temp1 = [[]]
    nr = 0
    for iline in range(len(line)):
        if line[iline] == "\n":
            temp1.append([])
            nr += 1
        else:
            temp1[nr].append(line[iline].split())
    #
    print("        Number of R : ", nr)
    norb = int(numpy.sqrt(len(temp1[0])-1)+0.1)
    print("    Number of bands : ", norb)

    irvec = numpy.zeros((nr, 3), numpy.int_)
    hopping = numpy.zeros((nr, norb, norb), numpy.complex_)

    for ir in range(nr):
        for ii in range(3):
            irvec[ir, ii] = int(temp1[ir][0][ii])

        ii = 0
        for iorb in range(norb):
            for jorb in range(norb):
                ii += 1
                hopping[ir, int(int(temp1[ir][ii][0]))-1, int(int(temp1[ir][ii][1]))-1] = \
                    float(temp1[ir][ii][2]) + 1.0j * float(temp1[ir][ii][3])
    #
    # Output to wannier90 format
    #
    with open(name_out, 'w') as f:

        print("Converted from RESPACK", file=f)
        print(norb, file=f)
        print(nr, file=f)
        for ir in range(nr):
            print("    1", end="", file=f)
            if ir % 15 == 14:
                print("", file=f)
        if nr % 15 != 0:
            print("", file=f)
        for ir in range(nr):
            for iorb in range(norb):
                for jorb in range(norb):
                    print("%5d%5d%5d%5d%5d%12.6f%12.6f" %
                          (irvec[ir, 0], irvec[ir, 1], irvec[ir, 2], jorb+1, iorb+1,
                           hopping[ir, jorb, iorb].real, hopping[ir, jorb, iorb].imag), file=f)


args = sys.argv

if len(args) != 2:
    print("\nUsage:\n")
    print("  $ respack2wan90.py seedname\n")
    exit(-1)

res2wan("./dir-wan/dat.h_mat_r", args[1] + "_hr.dat")
res2wan("./dir-intW/dat.Wmat", args[1] + "_ur.dat")
res2wan("./dir-intJ/dat.Jmat", args[1] + "_jr.dat")
