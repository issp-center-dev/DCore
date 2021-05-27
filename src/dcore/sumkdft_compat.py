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

from h5 import HDFArchive

def read_dft_input_data(file, subgrp, things_to_read):
    """

    Small version of SumkDFT.read_input_from_hdf()
    Read DFT data from a HDF file and return the data as a dict.

    """

    values = {}
    with HDFArchive(file, 'r') as ar:
        if not subgrp in ar:
            raise RuntimeError("subrp " + subgrp + "does not exist in " + file + "!")
        # first read the necessary things:
        for it in things_to_read:
            values[it] = ar[subgrp][it]

    return values


class SumkDFTCompat(object):
    """
    Reading data from a SumkDFT HDF5 file
    """
    def __init__(self, hdf_file, subgrp='dft_input'):

        things_to_read = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                          'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                          'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights',
                          'hopping',
                          'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']

        dft_data = read_dft_input_data(hdf_file, subgrp, things_to_read=things_to_read)

        for k, v in list(dft_data.items()):
            setattr(self, k, v)

        if self.SO != self.SP:
            raise RuntimeError("Not supported SP={} != SO={}.".format(self.SO, self.SP))