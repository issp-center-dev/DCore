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

from dcore.dcore_pre import dcore_pre
# from dcore._dispatcher import h5diff
from dcore.tools import h5diff
import os
#
# Execute dcore_pre.py to generate *test*.h5
# Then Check diff between *test*.h5 and the reference data *ref*.h5
#


def _uvals(interaction, l):
    if interaction == 'kanamori':
        return "kanamori = [(4.0, 2.2, 0.9)]"
    elif interaction == 'slater_uj':
        return f"slater_uj = [({l}, 4.0, 0.9)]"
    elif interaction == 'slater_f':
        return f"slater_f = [({l}, 4.0, 0.9, 0.1, 0.02)]"


# Sweep combinations of (interaction, l)
def test_interaction(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    for interaction in ['kanamori', 'slater_uj', 'slater_f']:
        for l in range(4):
            input_fname = f"interaction_{interaction}_{l}.in"
            seedname = f"interaction_{interaction}_{l}_test"
            seedname_ref = f"interaction_{interaction}_{l}_ref"

            with open(input_fname, 'w') as f:
                print("[model]", file=f)
                print(f"seedname = {seedname}", file=f)
                print(f"norb = {2*l+1}", file=f)
                print(f"interaction = {interaction}", file=f)
                print(_uvals(interaction, l), file=f)

            dcore_pre(input_fname)

            h5diff(seedname+".h5", seedname_ref+".h5", key='DCore/Umat')

    os.chdir(org_dir)


# Sweep combinations of (slater_basis, l)
def test_slater_basis(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    for slater_basis in ['spherical', 'spherical_j']:
        for l in range(4):
            input_fname = f"slater_basis_{slater_basis}_{l}.in"
            seedname = f"slater_basis_{slater_basis}_{l}_test"
            seedname_ref = f"slater_basis_{slater_basis}_{l}_ref"

            with open(input_fname, 'w') as f:
                print("[model]", file=f)
                print(f"seedname = {seedname}", file=f)
                print(f"norb = {2*l+1}", file=f)
                print(f"interaction = slater_uj", file=f)
                print(f"slater_uj = [({l}, 4.0, 0.9)]", file=f)
                print(f"slater_basis = {slater_basis}", file=f)

            dcore_pre(input_fname)

            h5diff(seedname+".h5", seedname_ref+".h5", key='DCore/Umat')

    os.chdir(org_dir)
