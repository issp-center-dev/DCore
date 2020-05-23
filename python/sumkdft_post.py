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

from .dft_tools_compat import SumkDFTTools
import pytriqs.utility.mpi as mpi

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


    # From dos_wannier_basis in dft_tools/python/sumk_dft_tools.py
    # This modified version allows to project DOS to any basis
    # Uses .data of only GfReFreq objects.
    def pdos(self, mu=None, broadening=None, mesh=None, with_Sigma=True, with_dc=True, save_to_file=True,
            rot_mat=None):
        """
        See the docstring of the original function in DFTTools

        Parameters
        ----------
        rot_mat : dict, optional
             Keys are strings representing block names.
             Values are 2D ndarrays representing the unitary matrix of a basis.
        """

        if (mesh is None) and (not with_Sigma):
            raise ValueError, "lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq."
        if mesh is None:
            om_mesh = [x.real for x in self.Sigma_imp_w[0].mesh]
            om_min = om_mesh[0]
            om_max = om_mesh[-1]
            n_om = len(om_mesh)
            mesh = (om_min, om_max, n_om)
        else:
            om_min, om_max, n_om = mesh
            om_mesh = numpy.linspace(om_min, om_max, n_om)

        G_loc = []
        for icrsh in range(self.n_corr_shells):
            spn = self.spin_block_names[self.corr_shells[icrsh]['SO']]
            glist = [GfReFreq(indices=inner, window=(om_min, om_max), n_points=n_om)
                     for block, inner in self.gf_struct_sumk[icrsh]]
            G_loc.append(
                BlockGf(name_list=spn, block_list=glist, make_copies=False))
        for icrsh in range(self.n_corr_shells):
            G_loc[icrsh].zero()

        DOS = {sp: numpy.zeros([n_om], numpy.float_)
               for sp in self.spin_block_names[self.SO]}
        DOSproj = [{} for ish in range(self.n_inequiv_shells)]
        DOSproj_orb = [{} for ish in range(self.n_inequiv_shells)]
        for ish in range(self.n_inequiv_shells):
            for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]:
                dim = self.corr_shells[self.inequiv_to_corr[ish]]['dim']
                DOSproj[ish][sp] = numpy.zeros([n_om], numpy.float_)
                DOSproj_orb[ish][sp] = numpy.zeros(
                    [n_om, dim, dim], numpy.complex_)

        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):

            G_latt_w = self.lattice_gf(
                ik=ik, mu=mu, iw_or_w="w", broadening=broadening, mesh=mesh, with_Sigma=with_Sigma, with_dc=with_dc)
            G_latt_w *= self.bz_weights[ik]

            # Non-projected DOS
            for iom in range(n_om):
                for bname, gf in G_latt_w:
                    DOS[bname][iom] -= gf.data[iom, :, :].imag.trace() / \
                        numpy.pi

            # Projected DOS:
            for icrsh in range(self.n_corr_shells):
                tmp = G_loc[icrsh].copy()
                for bname, gf in tmp:
                    tmp[bname] << self.downfold(ik, icrsh, bname, G_latt_w[
                                                bname], gf)  # downfolding G
                G_loc[icrsh] += tmp

        # Collect data from mpi:
        for bname in DOS:
            DOS[bname] = mpi.all_reduce(
                mpi.world, DOS[bname], lambda x, y: x + y)
        for icrsh in range(self.n_corr_shells):
            G_loc[icrsh] << mpi.all_reduce(
                mpi.world, G_loc[icrsh], lambda x, y: x + y)
        mpi.barrier()

        # Symmetrize and rotate to local coord. system if needed:
        if self.symm_op != 0:
            G_loc = self.symmcorr.symmetrize(G_loc)
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bname, gf in G_loc[icrsh]:
                    G_loc[icrsh][bname] << self.rotloc(
                        icrsh, gf, direction='toLocal')

        # ADDITION TO THE ORIGINAL FUNCTION
        #  Rotate basis
        if rot_mat is not None:
            for icrsh in range(self.n_corr_shells):
                for bname, gf in G_loc[icrsh]:
                    G_loc[icrsh][bname] << 
                        from_L_G_R(rot_mat[bname].transpose().conjugate(), G_loc[icrsh][bname], rot_mat[bname])

        # G_loc can now also be used to look at orbitally-resolved quantities
        for ish in range(self.n_inequiv_shells):
            for bname, gf in G_loc[self.inequiv_to_corr[ish]]:  # loop over spins
                for iom in range(n_om):
                    DOSproj[ish][bname][iom] -= gf.data[iom,
                                                        :, :].imag.trace() / numpy.pi
                DOSproj_orb[ish][bname][
                    :, :, :] += (1.0j*(gf-gf.conjugate().transpose())/2.0/numpy.pi).data[:,:,:]

        # Write to files
        if save_to_file and mpi.is_master_node():
            for sp in self.spin_block_names[self.SO]:
                f = open('DOS_wann_%s.dat' % sp, 'w')
                for iom in range(n_om):
                    f.write("%s    %s\n" % (om_mesh[iom], DOS[sp][iom]))
                f.close()

                # Partial
                for ish in range(self.n_inequiv_shells):
                    f = open('DOS_wann_%s_proj%s.dat' % (sp, ish), 'w')
                    for iom in range(n_om):
                        f.write("%s    %s\n" %
                                (om_mesh[iom], DOSproj[ish][sp][iom]))
                    f.close()

                    # Orbitally-resolved
                    for i in range(self.corr_shells[self.inequiv_to_corr[ish]]['dim']):
                        for j in range(i, self.corr_shells[self.inequiv_to_corr[ish]]['dim']):
                            f = open('DOS_wann_' + sp + '_proj' + str(ish) +
                                     '_' + str(i) + '_' + str(j) + '.dat', 'w')
                            for iom in range(n_om):
                                f.write("%s    %s    %s\n" % (
                                    om_mesh[iom], DOSproj_orb[ish][sp][iom, i, j].real,DOSproj_orb[ish][sp][iom, i, j].imag))
                            f.close()

        return DOS, DOSproj, DOSproj_orb
