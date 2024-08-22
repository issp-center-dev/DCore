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


def dcore_post(filename, np=1, prefix=None):
    """
    Removed function. Use dcore_anacont and dcore_spectrum instead.
    """

    sys.exit("dcore_post is removed. Use dcore_anacont and dcore_spectrum instead.")


def run():
    sys.exit("dcore_post is removed from DCore v4. Use dcore_anacont and dcore_spectrum instead.")
