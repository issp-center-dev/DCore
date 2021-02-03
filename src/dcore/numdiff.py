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

import sys


def numdiff(file1, file2, threshold=1.0e-6):
    """
    Main routine for the pre-processing tool

    Parameters
    ----------
    file1 : string
    file2 : string
        files to be compared
    """
    #
    # Input
    #
    print("Reading ", file1)
    dat1 = []
    with open(file1, 'r') as f:
        lines = f.readlines()
        for line in lines:
            nums = line.split()
            for num in nums:
                try:
                    dat1.append(float(num))
                except ValueError:
                    print("  String : ", num)
    print("Reading ", file2)
    dat2 = []
    with open(file2, 'r') as f:
        lines = f.readlines()
        for line in lines:
            nums = line.split()
            for num in nums:
                try:
                    dat2.append(float(num))
                except ValueError:
                    print("  String : ", num)

    print("Difference")
    if len(dat1) != len(dat2):
        print("DIFFERENT !  The length does not match.")
        exit(1)

    ave_diff = 0.0
    max_diff = 0.0
    for ii in range(len(dat1)):
        diff = abs(dat1[ii] - dat2[ii])
        ave_diff += diff
        if max_diff < diff:
            max_diff = diff

    ave_diff /= len(dat1)

    print("       Length = ", len(dat1))
    print("  Sum of diff = ", ave_diff)
    print("  Max of diff = ", max_diff)

    if max_diff < threshold:
        print("The difference is sufficiently small.")
    else:
        print("The difference is larger than the threshold.")
        exit(1)


if __name__ == '__main__':

    args = sys.argv

    if len(args) != 3:
        print("\nUsage:\n")
        print("  $ numdiff.py file1 file2\n")
        exit(-1)

    numdiff(args[1], args[2])
