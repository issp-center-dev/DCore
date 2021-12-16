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


import numpy

from triqs_dft_tools import SumkDFTTools
from .sumkdft_opt import SumkDFT_opt
import triqs.utility.mpi as mpi

class SumkDFTDCorePost(SumkDFTTools, SumkDFT_opt):
    """

    Extends the SumkDFTTools class with some tools for postprocessing in DCore.

    """

    def __init__(self, hdf_file, h_field=0.0, use_dft_blocks=False, dft_data='dft_input', symmcorr_data='dft_symmcorr_input',
                 parproj_data='dft_parproj_input', symmpar_data='dft_symmpar_input', bands_data='dft_bands_input',
                 transp_data='dft_transp_input', misc_data='dft_misc_input'):

        """
        Initialisation of the class. Parameters are exactly as for SumKDFT.
        """

        SumkDFTTools.__init__(self, hdf_file=hdf_file, h_field=h_field, use_dft_blocks=use_dft_blocks,
                         dft_data=dft_data, symmcorr_data=symmcorr_data, parproj_data=parproj_data,
                         symmpar_data=symmpar_data, bands_data=bands_data, transp_data=transp_data,
                         misc_data=misc_data)

    def calc_momentum_distribution(self, mu, beta, with_Sigma, with_dc):
        things_to_read = ['n_k', 'n_orbitals', 'proj_mat',
                          'hopping', 'n_parproj', 'proj_mat_all']
        value_read = self.read_input_from_hdf(
            subgrp=self.bands_data, things_to_read=things_to_read)

        den = numpy.zeros([self.n_k, 2-self.SO, self.n_orbitals[0, 0], self.n_orbitals[0, 0]], numpy.complex_)

        spn = self.spin_block_names[self.SO]

        ikarray = numpy.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):
            g_latt = self.lattice_gf(
                ik=ik, mu=mu, iw_or_w="iw", beta=beta, with_Sigma=with_Sigma, with_dc=with_dc)
            den0 = g_latt.density()
            g_latt.invert()
            for isp in range(len(spn)):
                den[ik, isp, :, :] = den0[spn[isp]][:, :]
                #
                # eigenvalue
                #
                #  ev[ik, isp, :] =  numpy.linalg.eigvalsh(g_latt[spn[isp]].data[len(g_latt[spn[isp]].mesh)//2, :, :])
        #
        # Collect density matrix across processes
        #
        den = mpi.all_reduce(mpi.world, den, lambda x, y: x + y)

        return den


# TODO: chi_0 etc.

