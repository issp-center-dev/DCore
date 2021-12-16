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

import argparse
import sys
import os

from triqs_cthyb import Solver as TRIQSCTHYBSolver
from .base import rotate_basis, make_h_int

from h5 import HDFArchive
import triqs.utility.mpi as mpi

from ..tools import convert_to_built_in_scalar_type


def main(input_file, output_file):
    """
    Solve the impurity problem.
    """

    import time
    t_start = time.time()

    with HDFArchive(os.path.abspath(input_file), 'r') as h:
        rot = h['rot'] if 'rot' in h else None
        beta = h['beta']
        gf_struct = h['gf_struct']
        # convert a dict to a list of pairs [ (str,[int,...]), ...]
        gf_struct = [(k, list(v)) for k, v in gf_struct.items()]
        u_mat = h['u_mat']
        n_iw = h['n_iw']
        G0_iw = h['G0_iw']
        params = h['params']

    use_spin_orbit = len(gf_struct) == 1

    # Rotate basis
    u_mat_rot = u_mat.copy()
    G0_iw_rot = G0_iw.copy()
    if not rot is None:
        u_mat_rot = rotate_basis(rot, use_spin_orbit, u_mat_rot, Gfs=[G0_iw_rot], direction='forward')
    h_int = make_h_int(u_mat_rot, gf_struct)

    # Create a working horse
    S = TRIQSCTHYBSolver(beta, gf_struct, n_iw)
    S.G0_iw << G0_iw_rot

    if 'random_seed' in params:
        raise RuntimeError("Do not set random_seed for triqs/cthyb manually")

    params['random_seed'] = 34788 + 928374 * mpi.rank + params['random_seed_offset']
    del params['random_seed_offset']

    for k in params:
        # e.g. numpy.bool_ to bool
        params[k] = convert_to_built_in_scalar_type(params[k])

    S.solve(h_int=h_int, **params)

    # Retrieve results from the working horse
    Sigma_iw = S.Sigma_iw.copy()
    G_iw = S.G_iw.copy()

    if not rot is None:
        Gfs = [Sigma_iw, G_iw]
        rotate_basis(rot, use_spin_orbit, u_matrix=None, Gfs=Gfs, direction='backward')

    if mpi.is_master_node():
        with HDFArchive(os.path.abspath(output_file), 'w') as h:
            h['Sigma_iw'] = Sigma_iw
            h['Gimp_iw'] = G_iw

    t_end = time.time()
    if mpi.is_master_node():
        print('TRIQS/cthyb ran for {} seconds.'.format(t_end-t_start))


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(description='Wrapper program for TRIQS/cthyb')
        parser.add_argument('input_file')
        parser.add_argument('output_file')
        args = parser.parse_args()
        main(args.input_file, args.output_file)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print("Unexpected error:", e)
        sys.exit(1)
