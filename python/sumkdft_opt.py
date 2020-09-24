from __future__ import print_function

import numpy
import copy
from pytriqs.applications.dft.sumk_dft import SumkDFT
from pytriqs.gf.local import *
import pytriqs.utility.mpi as mpi


class SumkDFT(SumkDFT):

    def __init__(self, hdf_file, h_field=0.0, use_dft_blocks=False,
                 dft_data='dft_input', symmcorr_data='dft_symmcorr_input', parproj_data='dft_parproj_input',
                 symmpar_data='dft_symmpar_input', bands_data='dft_bands_input', transp_data='dft_transp_input',
                 misc_data='dft_misc_input'):

        super(SumkDFT, self).__init__(hdf_file, h_field, use_dft_blocks, dft_data, symmcorr_data, parproj_data,
                                      symmpar_data, bands_data, transp_data, misc_data)

        # Store True if proj_mat is a unit matrix
        self.skip_projmat = self.check_if_proj_mat_is_a_unit_matrix()
        self.skip_projmat = False
        mpi.report("skip_projmat = {}".format(self.skip_projmat))

        # Make index for projection
        self.proj_index = self.make_proj_index()
        if self.proj_index is not None:
            mpi.report("proj_index has been set")

###############################################################
# OVERRIDE FUNCTIONS
# Modified parts are indicated by "+++MODIFIED"
# Added parts are indicated by "+++ADDED"
###############################################################

    def downfold(self, ik, ish, bname, gf_to_downfold, gf_inp, shells='corr', ir=None, overwrite_gf_inp=False,
                 fac=1.0):
        r"""
        If overwrite_gf_inp=True, the result is added to gf_inp, namely
            gf_inp += gf_downfolded * fac
        and None is returned.

        Parameters
        ----------
        overwrite_gf_inp : bool, optional

        fac : integer, optional

        """

        # +++ADDED
        if self.skip_projmat:
            if overwrite_gf_inp:
                gf_inp += gf_to_downfold * fac
                return None
            else:
                return gf_to_downfold.copy()
        if self.proj_index is not None:
            return self.downfold_index(ik, ish, bname, gf_to_downfold, gf_inp, shells, ir, overwrite_gf_inp, fac)

        gf_downfolded = gf_inp.copy()
        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
        elif shells == 'all':
            if ir is None:
                raise ValueError, "downfold: provide ir if treating all shells."
            dim = self.shells[ish]['dim']
            projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]

        gf_downfolded.from_L_G_R(
            projmat, gf_to_downfold, projmat.conjugate().transpose())

        # +++ADDED
        if overwrite_gf_inp:
            gf_inp += gf_downfolded * fac
            return None

        return gf_downfolded

    def upfold(self, ik, ish, bname, gf_to_upfold, gf_inp, shells='corr', ir=None, overwrite_gf_inp=False, fac=1.0):
        r"""
        If overwrite_gf_inp=True, the result is added to gf_inp, namely
            gf_inp += gf_upfolded * fac
        and None is returned.

        Parameters
        ----------
        overwrite_gf_inp : bool, optional

        fac : integer, optional

        """

        # +++ADDED
        if self.skip_projmat:
            if overwrite_gf_inp:
                gf_inp += gf_to_upfold * fac
                return None
            else:
                return gf_to_upfold.copy()
        if self.proj_index is not None:
            return self.upfold_index(ik, ish, bname, gf_to_upfold, gf_inp, shells, ir, overwrite_gf_inp, fac)

        gf_upfolded = gf_inp.copy()
        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
        elif shells == 'all':
            if ir is None:
                raise ValueError, "upfold: provide ir if treating all shells."
            dim = self.shells[ish]['dim']
            projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]

        gf_upfolded.from_L_G_R(
            projmat.conjugate().transpose(), gf_to_upfold, projmat)

        # +++ADDED
        if overwrite_gf_inp:
            gf_inp += gf_upfolded * fac
            return None

        return gf_upfolded

    def lattice_gf(self, ik, mu=None, iw_or_w="iw", beta=40, broadening=None, mesh=None, with_Sigma=True, with_dc=True):
        r"""
        """
        import time
        start = [time.time(),time.time()]
        def print_time(txt):
            # global start
            now = time.time()
            # print("time {:>8.5f}sec @ {}".format(now-start[0], txt))
            print("  time {:>8.5f}sec {:>8.5f}sec @ {}".format(now-start[0], now-start[1], txt))
            start[0] = now

        # print("\n=====================")
        # print("Start extract_G_loc")
        print_time("Start lattice_gf")


        if mu is None:
            mu = self.chemical_potential
        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        if (iw_or_w != "iw") and (iw_or_w != "w"):
            raise ValueError, "lattice_gf: Implemented only for Re/Im frequency functions."
        if not hasattr(self, "Sigma_imp_" + iw_or_w):
            with_Sigma = False
        if broadening is None:
            if mesh is None:
                broadening = 0.01
            else:  # broadening = 2 * \Delta omega, where \Delta omega is the spacing of omega points
                broadening = 2.0 * ((mesh[1] - mesh[0]) / (mesh[2] - 1))

        # Are we including Sigma?
        if with_Sigma:
            Sigma_imp = getattr(self, "Sigma_imp_" + iw_or_w)
            sigma_minus_dc = [s.copy() for s in Sigma_imp]
            if with_dc:
                sigma_minus_dc = self.add_dc(iw_or_w)
            if iw_or_w == "iw":
                # override beta if Sigma_iw is present
                beta = Sigma_imp[0].mesh.beta
                mesh = Sigma_imp[0].mesh
            elif iw_or_w == "w":
                mesh = Sigma_imp[0].mesh
                if broadening>0 and mpi.is_master_node():
                    warn('lattice_gf called with Sigma and broadening > 0 (broadening = {}). You might want to explicitly set the broadening to 0.'.format(broadening))
        else:
            if iw_or_w == "iw":
                if beta is None:
                    raise ValueError, "lattice_gf: Give the beta for the lattice GfReFreq."
                # Default number of Matsubara frequencies
                mesh = MeshImFreq(beta=beta, S='Fermion', n_max=1025)
            elif iw_or_w == "w":
                if mesh is None:
                    raise ValueError, "lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq."
                mesh = MeshReFreq(mesh[0], mesh[1], mesh[2])

        # Check if G_latt is present
        set_up_G_latt = False                       # Assume not
        if not hasattr(self, "G_latt_" + iw_or_w):
            # Need to create G_latt_(i)w
            set_up_G_latt = True
        else:                                       # Check that existing GF is consistent
            G_latt = getattr(self, "G_latt_" + iw_or_w)
            GFsize = [gf.N1 for bname, gf in G_latt]
            unchangedsize = all([self.n_orbitals[ik, ntoi[spn[isp]]] == GFsize[
                                isp] for isp in range(self.n_spin_blocks[self.SO])])
            if not unchangedsize:
                set_up_G_latt = True
            if (iw_or_w == "iw") and (self.G_latt_iw.mesh.beta != beta):
                set_up_G_latt = True  # additional check for ImFreq

        # Set up G_latt
        if set_up_G_latt:
            print("*** set_up_G_latt")
            block_structure = [
                range(self.n_orbitals[ik, ntoi[sp]]) for sp in spn]
            gf_struct = [(spn[isp], block_structure[isp])
                         for isp in range(self.n_spin_blocks[self.SO])]
            block_ind_list = [block for block, inner in gf_struct]
            if iw_or_w == "iw":
                glist = lambda: [GfImFreq(indices=inner, mesh=mesh)
                                 for block, inner in gf_struct]
            elif iw_or_w == "w":
                glist = lambda: [GfReFreq(indices=inner, mesh=mesh)
                                 for block, inner in gf_struct]
            G_latt = BlockGf(name_list=block_ind_list,
                             block_list=glist(), make_copies=False)
            G_latt.zero()

            # +++ADDED
            # Store iOmega_n or Omega to avoid re-computation
            self.cache_omega = BlockGf(name_list=block_ind_list,
                             block_list=glist(), make_copies=False)
            self.cache_omega.zero()
            if iw_or_w == "iw":
                self.cache_omega << iOmega_n
            elif iw_or_w == "w":
                self.cache_omega << Omega + 1j * broadening

        print_time("Set up G_latt")

        # if iw_or_w == "iw":
        #     G_latt << iOmega_n
        # elif iw_or_w == "w":
        #     G_latt << Omega + 1j * broadening
        # +++MODIFIED
        # just copy from cache
        G_latt << self.cache_omega

        print_time("G_latt << iOmega_n")

        idmat = [numpy.identity(
            self.n_orbitals[ik, ntoi[sp]], numpy.complex_) for sp in spn]
        M = copy.deepcopy(idmat)
        for ibl in range(self.n_spin_blocks[self.SO]):
            ind = ntoi[spn[ibl]]
            n_orb = self.n_orbitals[ik, ind]
            M[ibl] = self.hopping[ik, ind, 0:n_orb, 0:n_orb] - \
                (idmat[ibl] * mu) - (idmat[ibl] * self.h_field * (1 - 2 * ibl))
        G_latt -= M

        print_time("G_latt -= M")

        if with_Sigma:
            for icrsh in range(self.n_corr_shells):
                for bname, gf in G_latt:
                    # gf -= self.upfold(ik, icrsh, bname,
                    #                   sigma_minus_dc[icrsh][bname], gf)
                    # +++MODIFIED
                    # Subtract directly from G_latt (no temporary storage is introduced)
                    self.upfold(ik, icrsh, bname, sigma_minus_dc[icrsh][bname], gf, overwrite_gf_inp=True, fac=-1)

        print_time("upfold")

        G_latt.invert()
        setattr(self, "G_latt_" + iw_or_w, G_latt)
        print_time("invert")

        print_time("End lattice_gf")
        return G_latt

    def extract_G_loc(self, mu=None, iw_or_w='iw', with_Sigma=True, with_dc=True, broadening=None):
        r"""
        """

        import time
        start = [time.time(),time.time()]
        def print_time(txt):
            # global start
            now = time.time()
            # print("time {:>8.5f}sec @ {}".format(now-start[0], txt))
            print("time {:>8.5f}sec {:>8.5f}sec @ {}".format(now-start[0], now-start[1], txt))
            start[0] = now

        print("\n=====================")
        print("Start extract_G_loc")
        print_time("start")

        if mu is None:
            mu = self.chemical_potential

        if iw_or_w == "iw":
            G_loc = [self.Sigma_imp_iw[icrsh].copy() for icrsh in range(
                self.n_corr_shells)]   # this list will be returned
            beta = G_loc[0].mesh.beta
            G_loc_inequiv = [BlockGf(name_block_generator=[(block, GfImFreq(indices=inner, mesh=G_loc[0].mesh)) for block, inner in self.gf_struct_solver[ish].iteritems()],
                                     make_copies=False) for ish in range(self.n_inequiv_shells)]
        elif iw_or_w == "w":
            G_loc = [self.Sigma_imp_w[icrsh].copy() for icrsh in range(
                self.n_corr_shells)]   # this list will be returned
            mesh = G_loc[0].mesh
            G_loc_inequiv = [BlockGf(name_block_generator=[(block, GfReFreq(indices=inner, mesh=mesh)) for block, inner in self.gf_struct_solver[ish].iteritems()],
                                     make_copies=False) for ish in range(self.n_inequiv_shells)]

        for icrsh in range(self.n_corr_shells):
            G_loc[icrsh].zero()                          # initialize to zero

        print_time("k-sum start")

        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):
            print_time("in k-loop: k-sum")
            if iw_or_w == 'iw':
                G_latt = self.lattice_gf(
                    ik=ik, mu=mu, iw_or_w=iw_or_w, with_Sigma=with_Sigma, with_dc=with_dc, beta=beta)
            elif iw_or_w == 'w':
                mesh_parameters = (G_loc[0].mesh.omega_min,G_loc[0].mesh.omega_max,len(G_loc[0].mesh))
                G_latt = self.lattice_gf(
                    ik=ik, mu=mu, iw_or_w=iw_or_w, with_Sigma=with_Sigma, with_dc=with_dc, broadening=broadening, mesh=mesh_parameters)
            print_time("in k-loop: lattice_gf")
            G_latt *= self.bz_weights[ik]

            print_time("in k-loop: *bz_weights")
            for icrsh in range(self.n_corr_shells):
                # init temporary storage
                # tmp = G_loc[icrsh].copy()
                # for bname, gf in tmp:
                #     tmp[bname] << self.downfold(
                #         ik, icrsh, bname, G_latt[bname], gf)
                # G_loc[icrsh] += tmp
                # +++MODIFIED
                # Sum up directly into G_loc (no temporary storage is introduced)
                for bname, gf in G_loc[icrsh]:
                    self.downfold(ik, icrsh, bname, G_latt[bname], gf, overwrite_gf_inp=True, fac=+1)
            print_time("in k-loop: downfold")
        print_time("k-sum end")

        # Collect data from mpi
        for icrsh in range(self.n_corr_shells):
            G_loc[icrsh] << mpi.all_reduce(
                mpi.world, G_loc[icrsh], lambda x, y: x + y)
        mpi.barrier()
        print_time("mpi.all_reduce")

        # G_loc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:
        if self.symm_op != 0:
            G_loc = self.symmcorr.symmetrize(G_loc)

        # G_loc is rotated to the local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bname, gf in G_loc[icrsh]:
                    G_loc[icrsh][bname] << self.rotloc(
                        icrsh, gf, direction='toLocal')

        # transform to CTQMC blocks:
        for ish in range(self.n_inequiv_shells):
            for block, inner in self.gf_struct_solver[ish].iteritems():
                for ind1 in inner:
                    for ind2 in inner:
                        block_sumk, ind1_sumk = self.solver_to_sumk[
                            ish][(block, ind1)]
                        block_sumk, ind2_sumk = self.solver_to_sumk[
                            ish][(block, ind2)]
                        G_loc_inequiv[ish][block][ind1, ind2] << G_loc[
                            self.inequiv_to_corr[ish]][block_sumk][ind1_sumk, ind2_sumk]

        print_time("symm, rotations, solver_to_sumk")
        print("End extract_G_loc")
        print("=====================\n")

        # return only the inequivalent shells:
        return G_loc_inequiv

###############################################################
# ADDED FUNCTIONS
###############################################################

    def check_if_proj_mat_is_a_unit_matrix(self):
        """
        Return True if proj_mat is a unit matrix
        """

        def is_unitmatrix(mat):
            if mat.ndim == 2 \
                    and mat.shape[0] == mat.shape[1] \
                    and numpy.allclose(mat, numpy.identity(mat.shape[0])):
                return True
            else:
                return False

        # print(self.proj_mat.shape)
        projmat_is_a_unit_matrix = True
        for ik in range(self.n_k):
            for icrsh in range(self.n_corr_shells):
                for isp, sp in enumerate(self.spin_block_names[self.corr_shells[icrsh]['SO']]):
                    ind = self.spin_names_to_ind[
                        self.corr_shells[icrsh]['SO']][sp]
                    dim = self.corr_shells[icrsh]['dim']
                    n_orb = self.n_orbitals[ik, ind]
                    projmat = self.proj_mat[ik, ind, icrsh, 0:dim, 0:n_orb]
                    # print(projmat.shape)
                    # print(projmat)
                    if not is_unitmatrix(projmat):
                        projmat_is_a_unit_matrix = False
                        break

        return projmat_is_a_unit_matrix

    def make_proj_index(self):
        """
        Return (numpy.array)proj_index if proj_mat consists of only 0 and 1
               None otherwise
        """

        proj_index = numpy.zeros(self.proj_mat.shape[0:4], dtype=int)

        for ik in range(self.n_k):
            for icrsh in range(self.n_corr_shells):
                for isp, sp in enumerate(self.spin_block_names[self.corr_shells[icrsh]['SO']]):
                    ind = self.spin_names_to_ind[
                        self.corr_shells[icrsh]['SO']][sp]
                    dim = self.corr_shells[icrsh]['dim']
                    n_orb = self.n_orbitals[ik, ind]
                    projmat = self.proj_mat[ik, ind, icrsh, 0:dim, 0:n_orb]
                    # print(projmat.shape)
                    # print(projmat)
                    # index = numpy.where(projmat == 1)
                    # print(index)
                    for i in range(dim):
                    # for proj in projmat:
                        index = numpy.where(projmat[i] == 1)
                        # print(index)
                        if len(index) != 1:
                            return None
                        else:
                            proj_index[ik, ind, icrsh, i] = index[0]

        return proj_index

    def downfold_index(self, ik, ish, bname, gf_to_downfold, gf_inp, shells='corr', ir=None, overwrite_gf_inp=False, fac=1.0):
        r"""
        fancy-index version of downfold
        """

        # temp = gf_inp.copy()
        if overwrite_gf_inp:
            gf_downfolded = gf_inp
        else:
            gf_downfolded = gf_inp.copy()
            gf_downfolded.zero()
            gf_downfolded.tail.zero()

        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            # projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
            projindex = self.proj_index[ik, isp, ish, 0:dim]
        elif shells == 'all':
            if ir is None:
                raise ValueError, "downfold: provide ir if treating all shells."
            dim = self.shells[ish]['dim']
            # projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]
            projindex = self.proj_index_all[ik, isp, ish, ir, 0:dim]

        # gf_downfolded.from_L_G_R(
        #     projmat, gf_to_downfold, projmat.conjugate().transpose())
        # temp.from_L_G_R(
        #     projmat, gf_to_downfold, projmat.conjugate().transpose())

        # print("gf_downfolded.shape = {}".format(gf_downfolded.data.shape))
        # print("gf_to_downfold.shape = {}".format(gf_to_downfold.data.shape))

        # gf_downfolded.zero()
        i_start = projindex[0]
        i_end = projindex[0] + projindex.shape[0]
        # print(i_start, i_end)
        if numpy.allclose(projindex, list(range(i_start, i_end))):
            gf_downfolded.data[:, :, :] += gf_to_downfold.data[:, i_start:i_end, i_start:i_end] * fac
            gf_downfolded.tail.data[:, :, :] += gf_to_downfold.tail.data[:, i_start:i_end, i_start:i_end] * fac
        else:
            gf_downfolded.data[:, :, :] += gf_to_downfold.data[:, projindex, :][:, :, projindex] * fac
            gf_downfolded.tail.data[:, :, :] += gf_to_downfold.tail.data[:, projindex, :][:, :, projindex] * fac
        # print(gf_downfolded.tail)

        # assert numpy.allclose(gf_downfolded.data, temp.data)
        # assert numpy.allclose(gf_downfolded.tail.data, temp.tail.data)

        if overwrite_gf_inp:
            return None
        else:
            return gf_downfolded

    def upfold_index(self, ik, ish, bname, gf_to_upfold, gf_inp, shells='corr', ir=None, overwrite_gf_inp=False, fac=1.0):
        r"""
        fancy-index version of upfold
        """

        # temp = gf_inp.copy()
        if overwrite_gf_inp:
            gf_upfolded = gf_inp
        else:
            gf_upfolded = gf_inp.copy()
            gf_upfolded.zero()
            gf_upfolded.tail.zero()

        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            # projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
            projindex = self.proj_index[ik, isp, ish, 0:dim]
        elif shells == 'all':
            if ir is None:
                raise ValueError, "upfold: provide ir if treating all shells."
            dim = self.shells[ish]['dim']
            # projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]
            projindex = self.proj_index_all[ik, isp, ish, ir, 0:dim]

        # gf_upfolded.from_L_G_R(
        #     projmat.conjugate().transpose(), gf_to_upfold, projmat)
        # temp.from_L_G_R(
        #     projmat.conjugate().transpose(), gf_to_upfold, projmat)

        # print("gf_upfolded.shape = {}".format(gf_upfolded.data.shape))
        # print("gf_to_upfold.shape = {}".format(gf_to_upfold.data.shape))

        # print(projindex)
        i_start = projindex[0]
        i_end = projindex[0] + projindex.shape[0]
        # print(i_start, i_end)
        if numpy.allclose(projindex, list(range(i_start, i_end))):
        # if False:
            # print("continuous")
            # gf_upfolded.data[:, i_start:i_end, i_start:i_end] = gf_to_upfold.data
            # gf_upfolded.tail.data[:, i_start:i_end, i_start:i_end] = gf_to_upfold.tail.data
            # ***TEMP
            gf_upfolded.data[:, i_start:i_end, i_start:i_end] += gf_to_upfold.data * fac
            gf_upfolded.tail.data[:, i_start:i_end, i_start:i_end] += gf_to_upfold.tail.data * fac
        else:
            # gf_upfolded.data[:, projindex, :][:, :, projindex] = gf_to_upfold.data
            # gf_upfolded.tail.data[:, projindex, :][:, :, projindex] = gf_to_upfold.tail.data
            gf_upfolded.data[:, projindex, :][:, :, projindex] += gf_to_upfold.data * fac
            gf_upfolded.tail.data[:, projindex, :][:, :, projindex] += gf_to_upfold.tail.data * fac

        # projindex_inv = numpy.full(n_orb, dim, dtype=int)
        # for i, index in enumerate(projindex):
        #     projindex_inv[index] = i
        #
        # gf_to_upfold_data = numpy.zeros(gf_upfolded.data.shape, dtype=complex)
        # gf_to_upfold_data[:, 0:dim, 0:dim] = gf_to_upfold.data[:, :, :]
        # gf_upfolded.data[:, :, :] = gf_to_upfold_data[:, projindex_inv, :][:, :, projindex_inv]
        #
        # gf_to_upfold_tail = numpy.zeros(gf_upfolded.tail.data.shape, dtype=complex)
        # gf_to_upfold_tail[:, 0:dim, 0:dim] = gf_to_upfold.tail.data[:, :, :]
        # gf_upfolded.tail.data[:, :, :] = gf_to_upfold_tail[:, projindex_inv, :][:, :, projindex_inv]


        # gf_to_upfold_data_expanded = numpy.zeros(gf_upfolded.data.shape)
        # gf_to_upfold_data_expanded[:, 0:dim, 0:dim] = gf_to_upfold.data[:, :, :]


        # assert numpy.allclose(gf_upfolded.data, temp.data)
        # assert numpy.allclose(gf_upfolded.tail.data, temp.tail.data)

        if overwrite_gf_inp:
            return None
        else:
            return gf_upfolded
