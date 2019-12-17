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

import numpy

from .dft_tools_compat import SumkDFT, SumkDFTTools
import pytriqs.utility.mpi as mpi
import pytriqs.utility.dichotomy as dichotomy

class SumkDFTSCF(SumkDFT):
    """

    Change the behavior of calc_mu in SumkDFT

    """

    def __init__(self, hdf_file, h_field=0.0, use_dft_blocks=False, dft_data='dft_input', symmcorr_data='dft_symmcorr_input',
                 parproj_data='dft_parproj_input', symmpar_data='dft_symmpar_input', bands_data='dft_bands_input',
                 transp_data='dft_transp_input', misc_data='dft_misc_input'):

        """
        Initialisation of the class. Parameters are exactly as for SumKDFT.
        """

        SumkDFT.__init__(self, hdf_file=hdf_file, h_field=h_field, use_dft_blocks=use_dft_blocks,
                         dft_data=dft_data, symmcorr_data=symmcorr_data, parproj_data=parproj_data,
                         symmpar_data=symmpar_data, bands_data=bands_data, transp_data=transp_data,
                         misc_data=misc_data)


    def calc_mu(self, precision=0.01, iw_or_w='iw', broadening=None, delta=0.5, average=False):
        r"""
        Searches for the chemical potential that gives the DFT total charge.
        A simple bisection method is used.
        Parameters
        ----------
        precision : float, optional
                    A desired precision of the resulting total charge.
        iw_or_w : string, optional
                  - `iw_or_w` = 'iw' for a imaginary-frequency self-energy
                  - `iw_or_w` = 'w' for a real-frequency self-energy
        broadening : float, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
                     Only relevant for real-frequency GF.

        average : bool, optional
                     If True, compute chemical potential for target density +/- precision and average them.

        Returns
        -------
        mu : float
             Value of the chemical potential giving the DFT total charge
             within specified precision.
        """
        F = lambda mu: self.total_density(
            mu=mu, iw_or_w=iw_or_w, broadening=broadening).real
        density = self.density_required - self.charge_below

        if average:
            chem1 = dichotomy.dichotomy(function=F, x_init=self.chemical_potential, y_value=density-precision,
                                                          precision_on_y=precision, delta_x=delta, max_loops=100,
                                                          x_name="Chemical Potential", y_name="Total Density",
                                                          verbosity=3)[0]
            chem2 = dichotomy.dichotomy(function=F, x_init=chem1, y_value=density+precision,
                                        precision_on_y=precision, delta_x=delta, max_loops=100,
                                        x_name="Chemical Potential", y_name="Total Density",
                                        verbosity=3)[0]
            self.chemical_potential = 0.5*(chem1 + chem2)
        else:
            self.chemical_potential = dichotomy.dichotomy(function=F,
                                                      x_init=self.chemical_potential, y_value=density,
                                                      precision_on_y=precision, delta_x=delta, max_loops=100,
                                                      x_name="Chemical Potential", y_name="Total Density",
                                                      verbosity=3)[0]

        return self.chemical_potential



class SumkDFTDCorePost(SumkDFTTools):
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

        ikarray = numpy.array(range(self.n_k))
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
                #  ev[ik, isp, :] =  numpy.linalg.eigvalsh(g_latt[spn[isp]].data[len(g_latt[spn[isp]].mesh)/2, :, :])
        #
        # Collect density matrix across processes
        #
        den = mpi.all_reduce(mpi.world, den, lambda x, y: x + y)

        return den


# TODO: chi_0 etc.

