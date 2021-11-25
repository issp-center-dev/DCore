##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
from types import *
import numpy
from triqs.gf import *
from .. import mpi
from .symmetry import *
from .sumk_dft import SumkDFT
from scipy.integrate import *
from scipy.interpolate import *

if not hasattr(numpy, 'full'):
    # polyfill full for older numpy:
    numpy.full = lambda a, f: numpy.zeros(a) + f

class SumkDFTTools(SumkDFT):
    """
    Extends the SumkDFT class with some tools for analysing the data.
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

    # Uses .data of only GfReFreq objects.
    def dos_wannier_basis(self, mu=None, broadening=None, mesh=None, with_Sigma=True, with_dc=True, save_to_file=True):
        """
        Calculates the density of states in the basis of the Wannier functions.

        Parameters
        ----------
        mu : double, optional
             Chemical potential, overrides the one stored in the hdf5 archive.
        broadening : double, optional
                     Lorentzian broadening of the spectra. If not given, standard value of lattice_gf is used.
        mesh : real frequency MeshType, optional
               Omega mesh for the real-frequency Green's function. Given as parameter to lattice_gf.
        with_Sigma : boolean, optional
                     If True, the self energy is used for the calculation. If false, the DOS is calculated without self energy.
        with_dc : boolean, optional
                  If True the double counting correction is used.
        save_to_file : boolean, optional
                       If True, text files with the calculated data will be created.

        Returns
        -------
        DOS : Dict of numpy arrays
              Contains the full density of states.
        DOSproj :  Dict of numpy arrays
                   DOS projected to atoms.
        DOSproj_orb : Dict of numpy arrays
                      DOS projected to atoms and resolved into orbital contributions.
        """
        if (mesh is None) and (not with_Sigma):
            raise ValueError("lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq.")
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

        ikarray = numpy.array(list(range(self.n_k)))
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
    
    
    def dos_wannier_basis_all(self, mu=None, broadening=None, mesh=None, with_Sigma=True, with_dc=True, save_to_file=True):
        """
        Calculates the density of states in the basis of the Wannier functions.

        Parameters
        ----------
        mu : double, optional
             Chemical potential, overrides the one stored in the hdf5 archive.
        broadening : double, optional
                     Lorentzian broadening of the spectra. If not given, standard value of lattice_gf is used.
        mesh : real frequency MeshType, optional
               Omega mesh for the real-frequency Green's function. Given as parameter to lattice_gf.
        with_Sigma : boolean, optional
                     If True, the self energy is used for the calculation. If false, the DOS is calculated without self energy.
        with_dc : boolean, optional
                  If True the double counting correction is used.
        save_to_file : boolean, optional
                       If True, text files with the calculated data will be created.

        Returns
        -------
        DOS : Dict of numpy arrays
              Contains the full density of states.
        DOSproj :  Dict of numpy arrays
                   DOS projected to atoms.
        DOSproj_orb : Dict of numpy arrays
                      DOS projected to atoms and resolved into orbital contributions.
        """
        if (mesh is None) and (not with_Sigma):
            raise ValueError("lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq.")
        if mesh is None:
            om_mesh = [x.real for x in self.Sigma_imp_w[0].mesh]
            om_min = om_mesh[0]
            om_max = om_mesh[-1]
            n_om = len(om_mesh)
            mesh = (om_min, om_max, n_om)
        else:
            om_min, om_max, n_om = mesh
            om_mesh = numpy.linspace(om_min, om_max, n_om)

        spn = self.spin_block_names[self.SO]
        gf_struct_parproj = [[(sp, list(range(self.shells[ish]['dim']))) for sp in spn]
                             for ish in range(self.n_shells)]
        #print(self.proj_mat_csc.shape[2])
        #print(spn)
        n_local_orbs = self.proj_mat_csc.shape[2]
        gf_struct_parproj_all = [[(sp, list(range(n_local_orbs))) for sp in spn]]
    
        glist_all = [GfReFreq(indices=inner, window=(om_min, om_max), n_points=n_om)
                     for block, inner in gf_struct_parproj_all[0]]
        G_loc_all = BlockGf(name_list=spn, block_list=glist_all, make_copies=False)

        DOS = {sp: numpy.zeros([n_om], numpy.float_)
               for sp in self.spin_block_names[self.SO]}
        DOSproj = {}
        DOSproj_orb = {}

        for sp in self.spin_block_names[self.SO]:
            dim = n_local_orbs
            DOSproj[sp] = numpy.zeros([n_om], numpy.float_)
            DOSproj_orb[sp] = numpy.zeros(
                    [n_om, dim, dim], numpy.complex_)

        ikarray = numpy.array(list(range(self.n_k)))
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
            for bname, gf in G_latt_w:
                G_loc_all[bname] << self.downfold(ik, 0, bname, gf, G_loc_all[bname], shells='csc')
        # Collect data from mpi:
        for bname in DOS:
            DOS[bname] = mpi.all_reduce(
                mpi.world, DOS[bname], lambda x, y: x + y)
            G_loc_all[bname] << mpi.all_reduce(
                mpi.world, G_loc_all[bname], lambda x, y: x + y)
        mpi.barrier()

        # Symmetrize and rotate to local coord. system if needed:
        #if self.symm_op != 0:
        #    G_loc_all = self.symmcorr.symmetrize(G_loc_all)

        # G_loc can now also be used to look at orbitally-resolved quantities
        for bname, gf in G_loc_all:  # loop over spins
            for iom in range(n_om):
                DOSproj[bname][iom] -= gf.data[iom,:,:].imag.trace() / numpy.pi
            DOSproj_orb[bname][:,:,:] += (1.0j*(gf-gf.conjugate().transpose())/2.0/numpy.pi).data[:,:,:]
        # Write to files
        if save_to_file and mpi.is_master_node():
            for sp in self.spin_block_names[self.SO]:
                f = open('DOS_wann_%s.dat' % sp, 'w')
                for iom in range(n_om):
                    f.write("%s    %s\n" % (om_mesh[iom], DOS[sp][iom]))
                f.close()

                # Partial
                f = open('DOS_wann_all_%s_proj.dat' % (sp), 'w')
                for iom in range(n_om):
                    f.write("%s    %s\n" %
                            (om_mesh[iom], DOSproj[sp][iom]))
                f.close()

                # Orbitally-resolved
                for i in range(n_local_orbs):
                    for j in range(i, n_local_orbs):
                        f = open('DOS_wann_all' + sp + '_proj_' + str(i) + '_' + str(j) + '.dat', 'w')
                        for iom in range(n_om):
                                f.write("%s    %s    %s\n" % (
                                    om_mesh[iom], DOSproj_orb[sp][iom, i, j].real,DOSproj_orb[sp][iom, i, j].imag))
                        f.close()

        return DOS, DOSproj, DOSproj_orb

    # Uses .data of only GfReFreq objects.
    def dos_parproj_basis(self, mu=None, broadening=None, mesh=None, with_Sigma=True, with_dc=True, save_to_file=True):
        """
        Calculates the orbitally-resolved DOS.
        Different to dos_Wannier_basis is that here we calculate projections also to non-Wannier projectors, in the
        flavour of Wien2k QTL calculatuions.

        Parameters
        ----------
        mu : double, optional
             Chemical potential, overrides the one stored in the hdf5 archive.
        broadening : double, optional
                     Lorentzian broadening of the spectra. If not given, standard value of lattice_gf is used.
        mesh : real frequency MeshType, optional
               Omega mesh for the real-frequency Green's function. Given as parameter to lattice_gf.
        with_Sigma : boolean, optional
                     If True, the self energy is used for the calculation. If false, the DOS is calculated without self energy.
        with_dc : boolean, optional
                  If True the double counting correction is used.
        save_to_file : boolean, optional
                       If True, text files with the calculated data will be created.

        Returns
        -------
        DOS : Dict of numpy arrays
              Contains the full density of states.
        DOSproj :  Dict of numpy arrays
                   DOS projected to atoms.
        DOSproj_orb : Dict of numpy arrays
                      DOS projected to atoms and resolved into orbital contributions.
        """

        things_to_read = ['n_parproj', 'proj_mat_all',
                          'rot_mat_all', 'rot_mat_all_time_inv']
        value_read = self.read_input_from_hdf(
            subgrp=self.parproj_data, things_to_read=things_to_read)
        if not value_read:
            return value_read
        if self.symm_op:
            self.symmpar = Symmetry(self.hdf_file, subgroup=self.symmpar_data)

        if (mesh is None) and (not with_Sigma):
            raise ValueError("lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq.")
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
        spn = self.spin_block_names[self.SO]
        gf_struct_parproj = [[(sp, list(range(self.shells[ish]['dim']))) for sp in spn]
                             for ish in range(self.n_shells)]
        for ish in range(self.n_shells):
            glist = [GfReFreq(indices=inner, window=(om_min, om_max), n_points=n_om)
                     for block, inner in gf_struct_parproj[ish]]
            G_loc.append(
                BlockGf(name_list=spn, block_list=glist, make_copies=False))
        for ish in range(self.n_shells):
            G_loc[ish].zero()

        DOS = {sp: numpy.zeros([n_om], numpy.float_)
               for sp in self.spin_block_names[self.SO]}
        DOSproj = [{} for ish in range(self.n_shells)]
        DOSproj_orb = [{} for ish in range(self.n_shells)]
        for ish in range(self.n_shells):
            for sp in self.spin_block_names[self.SO]:
                dim = self.shells[ish]['dim']
                DOSproj[ish][sp] = numpy.zeros([n_om], numpy.float_)
                DOSproj_orb[ish][sp] = numpy.zeros(
                    [n_om, dim, dim], numpy.complex_)

        ikarray = numpy.array(list(range(self.n_k)))
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
            for ish in range(self.n_shells):
                tmp = G_loc[ish].copy()
                for ir in range(self.n_parproj[ish]):
                    for bname, gf in tmp:
                        tmp[bname] << self.downfold(ik, ish, bname, G_latt_w[
                                                    bname], gf, shells='all', ir=ir)
                    G_loc[ish] += tmp

        # Collect data from mpi:
        for bname in DOS:
            DOS[bname] = mpi.all_reduce(
                mpi.world, DOS[bname], lambda x, y: x + y)
        for ish in range(self.n_shells):
            G_loc[ish] << mpi.all_reduce(
                mpi.world, G_loc[ish], lambda x, y: x + y)
        mpi.barrier()

        # Symmetrize and rotate to local coord. system if needed:
        if self.symm_op != 0:
            G_loc = self.symmpar.symmetrize(G_loc)
        if self.use_rotations:
            for ish in range(self.n_shells):
                for bname, gf in G_loc[ish]:
                    G_loc[ish][bname] << self.rotloc(
                        ish, gf, direction='toLocal', shells='all')

        # G_loc can now also be used to look at orbitally-resolved quantities
        for ish in range(self.n_shells):
            for bname, gf in G_loc[ish]:
                for iom in range(n_om):
                    DOSproj[ish][bname][iom] -= gf.data[iom,
                                                        :, :].imag.trace() / numpy.pi
                DOSproj_orb[ish][bname][
                    :, :, :] += (1.0j*(gf-gf.conjugate().transpose())/2.0/numpy.pi).data[:,:,:]

        # Write to files
        if save_to_file and mpi.is_master_node():
            for sp in self.spin_block_names[self.SO]:
                f = open('DOS_parproj_%s.dat' % sp, 'w')
                for iom in range(n_om):
                    f.write("%s    %s\n" % (om_mesh[iom], DOS[sp][iom]))
                f.close()

                # Partial
                for ish in range(self.n_shells):
                    f = open('DOS_parproj_%s_proj%s.dat' % (sp, ish), 'w')
                    for iom in range(n_om):
                        f.write("%s    %s\n" %
                                (om_mesh[iom], DOSproj[ish][sp][iom]))
                    f.close()

                    # Orbitally-resolved
                    for i in range(self.shells[ish]['dim']):
                        for j in range(i, self.shells[ish]['dim']):
                            f = open('DOS_parproj_' + sp + '_proj' + str(ish) +
                                     '_' + str(i) + '_' + str(j) + '.dat', 'w')
                            for iom in range(n_om):
                                f.write("%s    %s    %s\n" % (
                                    om_mesh[iom], DOSproj_orb[ish][sp][iom, i, j].real,DOSproj_orb[ish][sp][iom, i, j].imag))
                            f.close()

        return DOS, DOSproj, DOSproj_orb

    # Uses .data of only GfReFreq objects.
    def spaghettis(self, broadening=None, plot_shift=0.0, plot_range=None, ishell=None, mu=None, save_to_file='Akw_'):
        """
        Calculates the correlated band structure using a real-frequency self energy.

        Parameters
        ----------
        mu : double, optional
             Chemical potential, overrides the one stored in the hdf5 archive.
        broadening : double, optional
                     Lorentzian broadening of the spectra. If not given, standard value of lattice_gf is used.
        plot_shift : double, optional
                     Offset for each A(k,w) for stacked plotting of spectra.
        plot_range : list of double, optional
                     Sets the energy window for plotting to (plot_range[0],plot_range[1]). If not provided, the energy mesh of the self energy is used.
        ishell : integer, optional
                 Contains the index of the shell on which the spectral function is projected. If ishell=None, the total spectrum without projection is calculated.
        save_to_file : string, optional
                       Filename where the spectra are stored.

        Returns
        -------
        Akw : Dict of numpy arrays
              Data as it is also written to the files.
        """

        assert hasattr(
            self, "Sigma_imp_w"), "spaghettis: Set Sigma_imp_w first."
        things_to_read = ['n_k', 'n_orbitals', 'proj_mat',
                          'hopping', 'n_parproj', 'proj_mat_all']
        value_read = self.read_input_from_hdf(
            subgrp=self.bands_data, things_to_read=things_to_read)
        if not value_read:
            return value_read
        if ishell is not None:
            things_to_read = ['rot_mat_all', 'rot_mat_all_time_inv']
            value_read = self.read_input_from_hdf(
                subgrp=self.parproj_data, things_to_read=things_to_read)
            if not value_read:
                return value_read

        if mu is None:
            mu = self.chemical_potential
        spn = self.spin_block_names[self.SO]
        mesh = [x.real for x in self.Sigma_imp_w[0].mesh]
        n_om = len(mesh)

        if plot_range is None:
            om_minplot = mesh[0] - 0.001
            om_maxplot = mesh[n_om - 1] + 0.001
        else:
            om_minplot = plot_range[0]
            om_maxplot = plot_range[1]

        if ishell is None:
            Akw = {sp: numpy.zeros([self.n_k, n_om], numpy.float_)
                   for sp in spn}
        else:
            Akw = {sp: numpy.zeros(
                [self.shells[ishell]['dim'], self.n_k, n_om], numpy.float_) for sp in spn}

        if not ishell is None:
            gf_struct_parproj = [
                (sp, list(range(self.shells[ishell]['dim']))) for sp in spn]
            G_loc = BlockGf(name_block_generator=[(block, GfReFreq(indices=inner, mesh=self.Sigma_imp_w[0].mesh))
                                                  for block, inner in gf_struct_parproj], make_copies=False)
            G_loc.zero()

        ikarray = numpy.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):

            G_latt_w = self.lattice_gf(
                ik=ik, mu=mu, iw_or_w="w", broadening=broadening)

            if ishell is None:
                # Non-projected A(k,w)
                for iom in range(n_om):
                    if (mesh[iom] > om_minplot) and (mesh[iom] < om_maxplot):
                        for bname, gf in G_latt_w:
                            Akw[bname][ik, iom] += gf.data[iom, :,
                                                           :].imag.trace() / (-1.0 * numpy.pi)
                        # shift Akw for plotting stacked k-resolved eps(k)
                        # curves
                        Akw[bname][ik, iom] += ik * plot_shift

            else:  # ishell not None
                # Projected A(k,w):
                G_loc.zero()
                tmp = G_loc.copy()
                for ir in range(self.n_parproj[ishell]):
                    for bname, gf in tmp:
                        tmp[bname] << self.downfold(ik, ishell, bname, G_latt_w[
                                                    bname], gf, shells='all', ir=ir)
                    G_loc += tmp

                # Rotate to local frame
                if self.use_rotations:
                    for bname, gf in G_loc:
                        G_loc[bname] << self.rotloc(
                            ishell, gf, direction='toLocal', shells='all')

                for iom in range(n_om):
                    if (mesh[iom] > om_minplot) and (mesh[iom] < om_maxplot):
                        for ish in range(self.shells[ishell]['dim']):
                            for sp in spn:
                                Akw[sp][ish, ik, iom] = G_loc[sp].data[
                                    iom, ish, ish].imag / (-1.0 * numpy.pi)

        # Collect data from mpi
        for sp in spn:
            Akw[sp] = mpi.all_reduce(mpi.world, Akw[sp], lambda x, y: x + y)
        mpi.barrier()

        if save_to_file and mpi.is_master_node():
            if ishell is None:
                for sp in spn:  # loop over GF blocs:
                    # Open file for storage:
                    f = open(save_to_file + sp + '.dat', 'w')
                    for ik in range(self.n_k):
                        for iom in range(n_om):
                            if (mesh[iom] > om_minplot) and (mesh[iom] < om_maxplot):
                                if plot_shift > 0.0001:
                                    f.write('%s      %s\n' %
                                            (mesh[iom], Akw[sp][ik, iom]))
                                else:
                                    f.write('%s     %s      %s\n' %
                                            (ik, mesh[iom], Akw[sp][ik, iom]))
                        f.write('\n')
                    f.close()

            else:  # ishell is not None
                for sp in spn:
                    for ish in range(self.shells[ishell]['dim']):
                        # Open file for storage:
                        f = open(save_to_file + str(ishell) + '_' +
                                 sp + '_proj' + str(ish) + '.dat', 'w')
                        for ik in range(self.n_k):
                            for iom in range(n_om):
                                if (mesh[iom] > om_minplot) and (mesh[iom] < om_maxplot):
                                    if plot_shift > 0.0001:
                                        f.write('%s      %s\n' % (
                                            mesh[iom], Akw[sp][ish, ik, iom]))
                                    else:
                                        f.write('%s     %s      %s\n' % (
                                            ik, mesh[iom], Akw[sp][ish, ik, iom]))
                            f.write('\n')
                        f.close()

        return Akw

    def partial_charges(self, beta=40, mu=None, with_Sigma=True, with_dc=True):
        """
        Calculates the orbitally-resolved density matrix for all the orbitals considered in the input, consistent with
        the definition of Wien2k. Hence, (possibly non-orthonormal) projectors have to be provided in the partial projectors subgroup of
        the hdf5 archive.

        Parameters
        ----------

        with_Sigma : boolean, optional
                     If True, the self energy is used for the calculation. If false, partial charges are calculated without self-energy correction.
        beta : double, optional
               In case the self-energy correction is not used, the inverse temperature where the calculation should be done has to be given here.
        mu : double, optional
             Chemical potential, overrides the one stored in the hdf5 archive.
        with_dc : boolean, optional
                  If True the double counting correction is used.

        Returns
        -------
        dens_mat : list of numpy array
                   A list of density matrices projected to all shells provided in the input.
        """

        things_to_read = ['dens_mat_below', 'n_parproj',
                          'proj_mat_all', 'rot_mat_all', 'rot_mat_all_time_inv']
        value_read = self.read_input_from_hdf(
            subgrp=self.parproj_data, things_to_read=things_to_read)
        if not value_read:
            return value_read
        if self.symm_op:
            self.symmpar = Symmetry(self.hdf_file, subgroup=self.symmpar_data)

        spn = self.spin_block_names[self.SO]
        ntoi = self.spin_names_to_ind[self.SO]
        # Density matrix in the window
        self.dens_mat_window = [[numpy.zeros([self.shells[ish]['dim'], self.shells[ish]['dim']], numpy.complex_)
                                 for ish in range(self.n_shells)]
                                for isp in range(len(spn))]
        # Set up G_loc
        gf_struct_parproj = [[(sp, list(range(self.shells[ish]['dim']))) for sp in spn]
                             for ish in range(self.n_shells)]
        if with_Sigma:
            G_loc = [BlockGf(name_block_generator=[(block, GfImFreq(indices=inner, mesh=self.Sigma_imp_iw[0].mesh))
                                                   for block, inner in gf_struct_parproj[ish]], make_copies=False)
                     for ish in range(self.n_shells)]
            beta = self.Sigma_imp_iw[0].mesh.beta
        else:
            G_loc = [BlockGf(name_block_generator=[(block, GfImFreq(indices=inner, beta=beta))
                                                   for block, inner in gf_struct_parproj[ish]], make_copies=False)
                     for ish in range(self.n_shells)]
        for ish in range(self.n_shells):
            G_loc[ish].zero()

        ikarray = numpy.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):

            G_latt_iw = self.lattice_gf(
                ik=ik, mu=mu, iw_or_w="iw", beta=beta, with_Sigma=with_Sigma, with_dc=with_dc)
            G_latt_iw *= self.bz_weights[ik]
            for ish in range(self.n_shells):
                tmp = G_loc[ish].copy()
                for ir in range(self.n_parproj[ish]):
                    for bname, gf in tmp:
                        tmp[bname] << self.downfold(ik, ish, bname, G_latt_iw[
                                                    bname], gf, shells='all', ir=ir)
                    G_loc[ish] += tmp

        # Collect data from mpi:
        for ish in range(self.n_shells):
            G_loc[ish] << mpi.all_reduce(
                mpi.world, G_loc[ish], lambda x, y: x + y)
        mpi.barrier()

        # Symmetrize and rotate to local coord. system if needed:
        if self.symm_op != 0:
            G_loc = self.symmpar.symmetrize(G_loc)
        if self.use_rotations:
            for ish in range(self.n_shells):
                for bname, gf in G_loc[ish]:
                    G_loc[ish][bname] << self.rotloc(
                        ish, gf, direction='toLocal', shells='all')

        for ish in range(self.n_shells):
            isp = 0
            for bname, gf in G_loc[ish]:
                self.dens_mat_window[isp][ish] = G_loc[ish].density()[bname]
                isp += 1

        # Add density matrices to get the total:
        dens_mat = [[self.dens_mat_below[ntoi[spn[isp]]][ish] + self.dens_mat_window[isp][ish]
                     for ish in range(self.n_shells)]
                    for isp in range(len(spn))]

        return dens_mat

    def print_hamiltonian(self):
        """
        Prints the Kohn-Sham Hamiltonian to the text files hamup.dat and hamdn.dat (no spin orbit-coupling), or to ham.dat (with spin-orbit coupling).
        """

        if self.SP == 1 and self.SO == 0:
            f1 = open('hamup.dat', 'w')
            f2 = open('hamdn.dat', 'w')
            for ik in range(self.n_k):
                for i in range(self.n_orbitals[ik, 0]):
                    f1.write('%s    %s\n' %
                             (ik, self.hopping[ik, 0, i, i].real))
                for i in range(self.n_orbitals[ik, 1]):
                    f2.write('%s    %s\n' %
                             (ik, self.hopping[ik, 1, i, i].real))
                f1.write('\n')
                f2.write('\n')
            f1.close()
            f2.close()
        else:
            f = open('ham.dat', 'w')
            for ik in range(self.n_k):
                for i in range(self.n_orbitals[ik, 0]):
                    f.write('%s    %s\n' %
                            (ik, self.hopping[ik, 0, i, i].real))
                f.write('\n')
            f.close()


# ----------------- transport -----------------------

    def read_transport_input_from_hdf(self):
        r"""
        Reads the data for transport calculations from the hdf5 archive.
        """
        thingstoread = ['band_window_optics', 'velocities_k']
        self.read_input_from_hdf(
            subgrp=self.transp_data, things_to_read=thingstoread)
        thingstoread = ['band_window', 'lattice_angles', 'lattice_constants',
                        'lattice_type', 'n_symmetries', 'rot_symmetries']
        self.read_input_from_hdf(
            subgrp=self.misc_data, things_to_read=thingstoread)

    def cellvolume(self, lattice_type, lattice_constants, latticeangle):
        r"""
        Determines the conventional und primitive unit cell volumes.

        Parameters
        ----------
        lattice_type : string
            Lattice type according to the Wien2k convention (P, F, B, R, H, CXY, CYZ, CXZ).
        lattice_constants : list of double
            Lattice constants (a, b, c).
        lattice angles : list of double
            Lattice angles (:math:`\alpha, \beta, \gamma`).

        Returns
        -------
        vol_c : double
            Conventional unit cell volume.
        vol_p : double
            Primitive unit cell volume.
        """

        a = lattice_constants[0]
        b = lattice_constants[1]
        c = lattice_constants[2]
        c_al = numpy.cos(latticeangle[0])
        c_be = numpy.cos(latticeangle[1])
        c_ga = numpy.cos(latticeangle[2])
        vol_c = a * b * c * \
            numpy.sqrt(1 + 2 * c_al * c_be * c_ga -
                       c_al ** 2 - c_be ** 2 - c_ga ** 2)

        det = {"P": 1, "F": 4, "B": 2, "R": 3,
               "H": 1, "CXY": 2, "CYZ": 2, "CXZ": 2}
        vol_p = vol_c / det[lattice_type]

        return vol_c, vol_p

    # Uses .data of only GfReFreq objects.
    def transport_distribution(self, beta, directions=['xx'], energy_window=None, Om_mesh=[0.0], with_Sigma=False, n_om=None, broadening=0.0):
        r"""
        Calculates the transport distribution 

        .. math::
           \Gamma_{\alpha\beta}\left(\omega+\Omega/2, \omega-\Omega/2\right) = \frac{1}{V} \sum_k Tr\left(v_{k,\alpha}A_{k}(\omega+\Omega/2)v_{k,\beta}A_{k}\left(\omega-\Omega/2\right)\right)

        in the direction :math:`\alpha\beta`. The velocities :math:`v_{k}` are read from the transport subgroup of the hdf5 archive. 

        Parameters
        ----------

        beta : double
            Inverse temperature :math:`\beta`.
        directions : list of double, optional
            :math:`\alpha\beta` e.g.: ['xx','yy','zz','xy','xz','yz'].
        energy_window : list of double, optional
            Specifies the upper and lower limit of the frequency integration for :math:`\Omega=0.0`. The window is automatically enlarged by the largest :math:`\Omega` value, 
            hence the integration is performed in the interval [energy_window[0]-max(Om_mesh), energy_window[1]+max(Om_mesh)].
        Om_mesh : list of double, optional
            :math:`\Omega` frequency mesh of the optical conductivity. For the conductivity and the Seebeck coefficient :math:`\Omega=0.0` has to be
            part of the mesh. In the current version Om_mesh is repined to the mesh provided by the self-energy! The actual mesh is printed on the screen and stored as
            member Om_mesh.
        with_Sigma : boolean, optional
            Determines whether the calculation is performed with or without self energy. If this parameter is set to False the self energy is set to zero (i.e. the DFT band 
            structure :math:`A(k,\omega)` is used). Note: For with_Sigma=False it is necessary to specify the parameters energy_window, n_om and broadening.
        n_om : integer, optional
            Number of equidistant frequency points in the interval [energy_window[0]-max(Om_mesh), energy_window[1]+max(Om_mesh)]. This parameters is only used if
            with_Sigma = False.
        broadening : double, optional
            Lorentzian broadening. It is necessary to specify the boradening if with_Sigma = False, otherwise this parameter can be set to 0.0.
        """

        # Check if wien converter was called and read transport subgroup form
        # hdf file
        if mpi.is_master_node():
            ar = HDFArchive(self.hdf_file, 'r')
            if not (self.transp_data in ar):
                raise IOError("transport_distribution: No %s subgroup in hdf file found! Call convert_transp_input first." % self.transp_data)
            # check if outputs file was converted
            if not ('n_symmetries' in ar['dft_misc_input']):
                raise IOError("transport_distribution: n_symmetries missing. Check if case.outputs file is present and call convert_misc_input() or convert_dft_input().")

        self.read_transport_input_from_hdf()

        if mpi.is_master_node():
            # k-dependent-projections.
            assert self.k_dep_projection == 1, "transport_distribution: k dependent projection is not implemented!"
            # positive Om_mesh
            assert all(
                Om >= 0.0 for Om in Om_mesh), "transport_distribution: Om_mesh should not contain negative values!"

        # Check if energy_window is sufficiently large and correct

        if (energy_window[0] >= energy_window[1] or energy_window[0] >= 0 or energy_window[1] <= 0):
            assert 0, "transport_distribution: energy_window wrong!"

        if (abs(self.fermi_dis(energy_window[0], beta) * self.fermi_dis(-energy_window[0], beta)) > 1e-5
                or abs(self.fermi_dis(energy_window[1], beta) * self.fermi_dis(-energy_window[1], beta)) > 1e-5):
            mpi.report(
                "\n####################################################################")
            mpi.report(
                "transport_distribution: WARNING - energy window might be too narrow!")
            mpi.report(
                "####################################################################\n")

        # up and down are equivalent if SP = 0
        n_inequiv_spin_blocks = self.SP + 1 - self.SO
        self.directions = directions
        dir_to_int = {'x': 0, 'y': 1, 'z': 2}

        # calculate A(k,w)
        #######################################

        # Define mesh for Green's function and in the specified energy window
        if (with_Sigma == True):
            self.omega = numpy.array([round(x.real, 12)
                                      for x in self.Sigma_imp_w[0].mesh])
            mesh = None
            mu = self.chemical_potential
            n_om = len(self.omega)
            mpi.report("Using omega mesh provided by Sigma!")

            if energy_window is not None:
                # Find according window in Sigma mesh
                ioffset = numpy.sum(
                    self.omega < energy_window[0] - max(Om_mesh))
                self.omega = self.omega[numpy.logical_and(self.omega >= energy_window[
                                                          0] - max(Om_mesh), self.omega <= energy_window[1] + max(Om_mesh))]
                n_om = len(self.omega)

                # Truncate Sigma to given omega window
                # In the future there should be an option in gf to manipulate the mesh (e.g. truncate) directly.
                # For now we stick with this:
                for icrsh in range(self.n_corr_shells):
                    Sigma_save = self.Sigma_imp_w[icrsh].copy()
                    spn = self.spin_block_names[self.corr_shells[icrsh]['SO']]
                    glist = lambda: [GfReFreq(indices=inner, window=(self.omega[
                                              0], self.omega[-1]), n_points=n_om) for block, inner in self.gf_struct_sumk[icrsh]]
                    self.Sigma_imp_w[icrsh] = BlockGf(
                        name_list=spn, block_list=glist(), make_copies=False)
                    for i, g in self.Sigma_imp_w[icrsh]:
                        for iL in g.indices[0]:
                            for iR in g.indices[0]:
                                for iom in range(n_om):
                                    g.data[iom, int(iL), int(iR)] = Sigma_save[
                                        i].data[ioffset + iom, int(iL), int(iR)]
        else:
            assert n_om is not None, "transport_distribution: Number of omega points (n_om) needed to calculate transport distribution!"
            assert energy_window is not None, "transport_distribution: Energy window needed to calculate transport distribution!"
            assert broadening != 0.0 and broadening is not None, "transport_distribution: Broadening necessary to calculate transport distribution!"
            self.omega = numpy.linspace(
                energy_window[0] - max(Om_mesh), energy_window[1] + max(Om_mesh), n_om)
            mesh = [energy_window[0] -
                    max(Om_mesh), energy_window[1] + max(Om_mesh), n_om]
            mu = 0.0

        # Define mesh for optic conductivity
        d_omega = round(numpy.abs(self.omega[0] - self.omega[1]), 12)
        iOm_mesh = numpy.array([round((Om / d_omega), 0) for Om in Om_mesh])
        self.Om_mesh = iOm_mesh * d_omega

        if mpi.is_master_node():
            print("Chemical potential: ", mu)
            print("Using n_om = %s points in the energy_window [%s,%s]" % (n_om, self.omega[0], self.omega[-1]), end=' ')
            print("where the omega vector is:")
            print(self.omega)
            print("Calculation requested for Omega mesh:   ", numpy.array(Om_mesh))
            print("Omega mesh automatically repined to:  ", self.Om_mesh)

        self.Gamma_w = {direction: numpy.zeros(
            (len(self.Om_mesh), n_om), dtype=numpy.float_) for direction in self.directions}

        # Sum over all k-points
        ikarray = numpy.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):
            # Calculate G_w  for ik and initialize A_kw
            G_w = self.lattice_gf(ik, mu, iw_or_w="w", beta=beta,
                                  broadening=broadening, mesh=mesh, with_Sigma=with_Sigma)
            A_kw = [numpy.zeros((self.n_orbitals[ik][isp], self.n_orbitals[ik][isp], n_om), dtype=numpy.complex_)
                    for isp in range(n_inequiv_spin_blocks)]

            for isp in range(n_inequiv_spin_blocks):
                # copy data from G_w (swapaxes is used to have omega in the 3rd
                # dimension)
                A_kw[isp] = copy.deepcopy(G_w[self.spin_block_names[self.SO][
                                          isp]].data.swapaxes(0, 1).swapaxes(1, 2))
                # calculate A(k,w) for each frequency
                for iw in range(n_om):
                    A_kw[isp][:, :, iw] = -1.0 / (2.0 * numpy.pi * 1j) * (
                        A_kw[isp][:, :, iw] - numpy.conjugate(numpy.transpose(A_kw[isp][:, :, iw])))

                b_min = max(self.band_window[isp][
                            ik, 0], self.band_window_optics[isp][ik, 0])
                b_max = min(self.band_window[isp][
                            ik, 1], self.band_window_optics[isp][ik, 1])
                A_i = slice(
                    b_min - self.band_window[isp][ik, 0], b_max - self.band_window[isp][ik, 0] + 1)
                v_i = slice(b_min - self.band_window_optics[isp][
                            ik, 0], b_max - self.band_window_optics[isp][ik, 0] + 1)

                # loop over all symmetries
                for R in self.rot_symmetries:
                    # get transformed velocity under symmetry R
                    vel_R = copy.deepcopy(self.velocities_k[isp][ik])
                    for nu1 in range(self.band_window_optics[isp][ik, 1] - self.band_window_optics[isp][ik, 0] + 1):
                        for nu2 in range(self.band_window_optics[isp][ik, 1] - self.band_window_optics[isp][ik, 0] + 1):
                            vel_R[nu1][nu2][:] = numpy.dot(
                                R, vel_R[nu1][nu2][:])

                    # calculate Gamma_w for each direction from the velocities
                    # vel_R and the spectral function A_kw
                    for direction in self.directions:
                        for iw in range(n_om):
                            for iq in range(len(self.Om_mesh)):
                                if(iw + iOm_mesh[iq] >= n_om or self.omega[iw] < -self.Om_mesh[iq] + energy_window[0] or self.omega[iw] > self.Om_mesh[iq] + energy_window[1]):
                                    continue

                                self.Gamma_w[direction][iq, iw] += (numpy.dot(numpy.dot(numpy.dot(vel_R[v_i, v_i, dir_to_int[direction[0]]],
                                                                                                  A_kw[isp][A_i, A_i, int(iw + iOm_mesh[iq])]), vel_R[v_i, v_i, dir_to_int[direction[1]]]),
                                                                              A_kw[isp][A_i, A_i, iw]).trace().real * self.bz_weights[ik])

        for direction in self.directions:
            self.Gamma_w[direction] = (mpi.all_reduce(mpi.world, self.Gamma_w[direction], lambda x, y: x + y)
                                       / self.cellvolume(self.lattice_type, self.lattice_constants, self.lattice_angles)[1] / self.n_symmetries)

    def transport_coefficient(self, direction, iq, n, beta, method=None):
        r"""
        Calculates the transport coefficient A_n in a given direction for a given :math:`\Omega`. The required members (Gamma_w, directions, Om_mesh) have to be obtained first
        by calling the function :meth:`transport_distribution <dft.sumk_dft_tools.SumkDFTTools.transport_distribution>`. For n>0 A is set to NaN if :math:`\Omega` is not 0.0. 

        Parameters
        ----------
        direction : string
           :math:`\alpha\beta` e.g.: 'xx','yy','zz','xy','xz','yz'.
        iq : integer
            Index of :math:`\Omega` point in the member Om_mesh.
        n : integer
            Number of the desired moment of the transport distribution.
        beta : double
            Inverse temperature :math:`\beta`.
        method : string
            Integration method: cubic spline and scipy.integrate.quad ('quad'), simpson rule ('simps'), trapezoidal rule ('trapz'), rectangular integration (otherwise)
            Note that the sampling points of the the self-energy are used!

        Returns
        -------
        A : double
            Transport coefficient.
        """

        if not (mpi.is_master_node()):
            return

        assert hasattr(
            self, 'Gamma_w'), "transport_coefficient: Run transport_distribution first or load data from h5!"

        if (self.Om_mesh[iq] == 0.0 or n == 0.0):
            A = 0.0
            # setup the integrand
            if (self.Om_mesh[iq] == 0.0):
                A_int = self.Gamma_w[direction][iq] * (self.fermi_dis(
                    self.omega, beta) * self.fermi_dis(-self.omega, beta)) * (self.omega * beta)**n
            elif (n == 0.0):
                A_int = self.Gamma_w[direction][iq] * (self.fermi_dis(self.omega, beta) - self.fermi_dis(
                    self.omega + self.Om_mesh[iq], beta)) / (self.Om_mesh[iq] * beta)

            # w-integration
            if method == 'quad':
                # quad on interpolated w-points with cubic spline
                A_int_interp = interp1d(self.omega, A_int, kind='cubic')
                A = quad(A_int_interp, min(self.omega), max(self.omega),
                         epsabs=1.0e-12, epsrel=1.0e-12, limit=500)
                A = A[0]
            elif method == 'simps':
                # simpson rule for w-grid
                A = simps(A_int, self.omega)
            elif method == 'trapz':
                # trapezoidal rule for w-grid
                A = numpy.trapz(A_int, self.omega)
            else:
                # rectangular integration for w-grid (orignal implementation)
                d_w = self.omega[1] - self.omega[0]
                for iw in range(self.Gamma_w[direction].shape[1]):
                    A += A_int[iw] * d_w
            A = A * numpy.pi * (2.0 - self.SP)
        else:
            A = numpy.nan
        return A

    def conductivity_and_seebeck(self, beta, method=None):
        r"""
        Calculates the Seebeck coefficient and the optical conductivity by calling 
        :meth:`transport_coefficient <dft.sumk_dft_tools.SumkDFTTools.transport_coefficient>`. 
        The required members (Gamma_w, directions, Om_mesh) have to be obtained first by calling the function 
        :meth:`transport_distribution <dft.sumk_dft_tools.SumkDFTTools.transport_distribution>`. 

        Parameters
        ----------
        beta : double
            Inverse temperature :math:`\beta`.

        Returns
        -------
        optic_cond : dictionary of double vectors
            Optical conductivity in each direction and frequency given by Om_mesh.

        seebeck : dictionary of double
            Seebeck coefficient in each direction. If zero is not present in Om_mesh the Seebeck coefficient is set to NaN.

        kappa : dictionary of double. 
            thermal conductivity in each direction. If zero is not present in Om_mesh the thermal conductivity is set to NaN
        """

        if not (mpi.is_master_node()):
            return

        assert hasattr(
            self, 'Gamma_w'), "conductivity_and_seebeck: Run transport_distribution first or load data from h5!"
        n_q = self.Gamma_w[self.directions[0]].shape[0]

        A0 = {direction: numpy.full((n_q,), numpy.nan)
              for direction in self.directions}
        A1 = {direction: numpy.full((n_q,), numpy.nan)
              for direction in self.directions}
        A2 = {direction: numpy.full((n_q,), numpy.nan)
              for direction in self.directions}
        self.seebeck = {direction: numpy.nan for direction in self.directions}
        self.kappa = {direction: numpy.nan for direction in self.directions}
        self.optic_cond = {direction: numpy.full(
            (n_q,), numpy.nan) for direction in self.directions}

        for direction in self.directions:
            for iq in range(n_q):
                A0[direction][iq] = self.transport_coefficient(
                    direction, iq=iq, n=0, beta=beta, method=method)
                A1[direction][iq] = self.transport_coefficient(
                    direction, iq=iq, n=1, beta=beta, method=method)
                A2[direction][iq] = self.transport_coefficient(
                    direction, iq=iq, n=2, beta=beta, method=method)
                print("A_0 in direction %s for Omega = %.2f    %e a.u." % (direction, self.Om_mesh[iq], A0[direction][iq]))
                print("A_1 in direction %s for Omega = %.2f    %e a.u." % (direction, self.Om_mesh[iq], A1[direction][iq]))
                print("A_2 in direction %s for Omega = %.2f    %e a.u." % (direction, self.Om_mesh[iq], A2[direction][iq]))
                if ~numpy.isnan(A1[direction][iq]):
                    # Seebeck and kappa are overwritten if there is more than one Omega =
                    # 0 in Om_mesh
                    self.seebeck[direction] = - \
                        A1[direction][iq] / A0[direction][iq] * 86.17
                    self.kappa[direction] = A2[direction][iq] - A1[direction][iq]*A1[direction][iq]/A0[direction][iq]
                    self.kappa[direction] *= 293178.0
            self.optic_cond[direction] = beta * \
                A0[direction] * 10700.0 / numpy.pi
            for iq in range(n_q):
                print("Conductivity in direction %s for Omega = %.2f       %f  x 10^4 Ohm^-1 cm^-1" % (direction, self.Om_mesh[iq], self.optic_cond[direction][iq]))
                if not (numpy.isnan(A1[direction][iq])):
                    print("Seebeck in direction      %s for Omega = 0.00      %f  x 10^(-6) V/K" % (direction, self.seebeck[direction]))
                    print("kappa in direction      %s for Omega = 0.00      %f  W/(m * K)" % (direction, self.kappa[direction]))

        return self.optic_cond, self.seebeck, self.kappa


    def fermi_dis(self, w, beta):
        r"""
        Fermi distribution.

        .. math::
           f(x) = 1/(e^x+1).

        Parameters
        ----------
        w : double
           frequency
        beta : double
           inverse temperature

        Returns
        -------
        f : double
        """
        return 1.0 / (numpy.exp(w * beta) + 1)
