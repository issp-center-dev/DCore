import numpy
import copy
from warnings import warn
from dcore._dispatcher import *
from .tools import calc_total_density, calc_density_matrix

class SumkDFT_opt(SumkDFT):

    def __init__(self, hdf_file, h_field=0.0, use_dft_blocks=False,
                 dft_data='dft_input', symmcorr_data='dft_symmcorr_input', parproj_data='dft_parproj_input',
                 symmpar_data='dft_symmpar_input', bands_data='dft_bands_input', transp_data='dft_transp_input',
                 misc_data='dft_misc_input'):

        super().__init__(hdf_file, h_field, use_dft_blocks, dft_data, symmcorr_data, parproj_data,
                                      symmpar_data, bands_data, transp_data, misc_data)

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

        fac : float, optional

        """

        # +++ADDED
        if self.index_works(shells):
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
                raise ValueError("downfold: provide ir if treating all shells.")
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

        fac : float, optional

        """

        # +++ADDED
        if self.index_works(shells):
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
                raise ValueError("upfold: provide ir if treating all shells.")
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

        # +++ADDED
        # print_time function is inserted in several places for benchmark
        import time
        start = [time.time(),time.time()]
        def print_time(txt):
            # now = time.time()
            # print("  time {:>8.5f}sec {:>8.5f}sec @ {}".format(now-start[0], now-start[1], txt))
            # start[0] = now
            pass
        print_time("Start lattice_gf")

        if mu is None:
            mu = self.chemical_potential
        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        if (iw_or_w != "iw") and (iw_or_w != "w"):
            raise ValueError("lattice_gf: Implemented only for Re/Im frequency functions.")
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
                    raise ValueError("lattice_gf: Give the beta for the lattice GfReFreq.")
                # Default number of Matsubara frequencies
                mesh = MeshImFreq(beta=beta, S='Fermion', n_max=1025)
            elif iw_or_w == "w":
                if mesh is None:
                    raise ValueError("lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq.")
                mesh = MeshReFreq(mesh[0], mesh[1], mesh[2])

        # Check if G_latt is present
        set_up_G_latt = False                       # Assume not
        if not hasattr(self, "G_latt_" + iw_or_w):
            # Need to create G_latt_(i)w
            set_up_G_latt = True
        else:                                       # Check that existing GF is consistent
            G_latt = getattr(self, "G_latt_" + iw_or_w)
            GFsize = [gf.data.shape[1] for bname, gf in G_latt]
            unchangedsize = all([self.n_orbitals[ik, ntoi[spn[isp]]] == GFsize[
                                isp] for isp in range(self.n_spin_blocks[self.SO])])
            if not unchangedsize:
                set_up_G_latt = True
            if (iw_or_w == "iw") and (self.G_latt_iw.mesh.beta != beta):
                set_up_G_latt = True  # additional check for ImFreq

        # Set up G_latt
        if set_up_G_latt:
            block_structure = [
                list(range(self.n_orbitals[ik, ntoi[sp]])) for sp in spn]
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

        # +++ADDED
        # print_time function is inserted in several places for benchmark
        import time
        start = [time.time(),time.time()]
        def print_time(txt):
            # now = time.time()
            # print("time {:>8.5f}sec {:>8.5f}sec @ {}".format(now-start[0], now-start[1], txt))
            # start[0] = now
            pass
        # print("\n=====================")
        # print("Start extract_G_loc")
        print_time("start")

        if mu is None:
            mu = self.chemical_potential

        if iw_or_w == "iw":
            G_loc = [self.Sigma_imp_iw[icrsh].copy() for icrsh in range(
                self.n_corr_shells)]   # this list will be returned
            beta = G_loc[0].mesh.beta
            G_loc_inequiv = [BlockGf(name_block_generator=[(block, GfImFreq(indices=inner, mesh=G_loc[0].mesh)) for block, inner in self.gf_struct_solver[ish].items()],
                                     make_copies=False) for ish in range(self.n_inequiv_shells)]
        elif iw_or_w == "w":
            G_loc = [self.Sigma_imp_w[icrsh].copy() for icrsh in range(
                self.n_corr_shells)]   # this list will be returned
            mesh = G_loc[0].mesh
            G_loc_inequiv = [BlockGf(name_block_generator=[(block, GfReFreq(indices=inner, mesh=mesh)) for block, inner in self.gf_struct_solver[ish].items()],
                                     make_copies=False) for ish in range(self.n_inequiv_shells)]

        for icrsh in range(self.n_corr_shells):
            G_loc[icrsh].zero()                          # initialize to zero

        print_time("k-sum start")

        ikarray = numpy.array(list(range(self.n_k)))
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
            for block, inner in self.gf_struct_solver[ish].items():
                for ind1 in inner:
                    for ind2 in inner:
                        block_sumk, ind1_sumk = self.solver_to_sumk[
                            ish][(block, ind1)]
                        block_sumk, ind2_sumk = self.solver_to_sumk[
                            ish][(block, ind2)]
                        G_loc_inequiv[ish][block][ind1, ind2] << G_loc[
                            self.inequiv_to_corr[ish]][block_sumk][ind1_sumk, ind2_sumk]

        print_time("symm, rotations, solver_to_sumk")
        # print("End extract_G_loc")
        # print("=====================\n")

        # return only the inequivalent shells:
        return G_loc_inequiv

    ###############################################################
    # ADDED FUNCTIONS
    ###############################################################

    def index_works(self, shells):
        """
        Return True if the index version of upfold/downfold can be used.
        If True, proj_index is set (necessary for running the index version).
        """
        if shells=='corr':
            # (Re)make proj_index if proj_index has not been computed or if proj_mat has been updated
            if not hasattr(self, 'projmat_id') or self.projmat_id != id(self.proj_mat):
                self.proj_index = self.make_proj_index()
                self.projmat_id = id(self.proj_mat)  # Store ID(proj_mat) to skip making the index next time
                if self.proj_index is not None:
                    mpi.report("The fancy-index version of upfold/downfold is used.")
            projindex = self.proj_index
        elif shells=='all':
            # TODO: make self.proj_index_all from self.proj_mat_all
            projindex = None

        return projindex is not None

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
                    for i in range(dim):
                        index = numpy.where(projmat[i] == 1)
                        if len(index) != 1:
                            return None
                        else:
                            proj_index[ik, ind, icrsh, i] = index[0]

        return proj_index

    def downfold_index(self, ik, ish, bname, gf_to_downfold, gf_inp, shells='corr', ir=None, overwrite_gf_inp=False, fac=1.0):
        r"""
        fancy-index version of downfold
        """

        if overwrite_gf_inp:
            gf_downfolded = gf_inp
        else:
            gf_downfolded = gf_inp.copy()
            gf_downfolded.zero()

        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            # projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
            projindex = self.proj_index[ik, isp, ish, 0:dim]
        elif shells == 'all':
            if ir is None:
                raise ValueError("downfold: provide ir if treating all shells.")
            dim = self.shells[ish]['dim']
            # projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]
            projindex = self.proj_index_all[ik, isp, ish, ir, 0:dim]

        # gf_downfolded.from_L_G_R(
        #     projmat, gf_to_downfold, projmat.conjugate().transpose())

        i_start = projindex[0]
        i_end = projindex[0] + projindex.shape[0]
        if numpy.allclose(projindex, list(range(i_start, i_end))):
            gf_downfolded.data[:, :, :] += gf_to_downfold.data[:, i_start:i_end, i_start:i_end] * fac
        else:
            gf_downfolded.data[:, :, :] += gf_to_downfold.data[:, projindex, :][:, :, projindex] * fac

        if overwrite_gf_inp:
            return None
        else:
            return gf_downfolded

    def upfold_index(self, ik, ish, bname, gf_to_upfold, gf_inp, shells='corr', ir=None, overwrite_gf_inp=False, fac=1.0):
        r"""
        fancy-index version of upfold
        """

        if overwrite_gf_inp:
            gf_upfolded = gf_inp
        else:
            gf_upfolded = gf_inp.copy()
            gf_upfolded.zero()

        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            # projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
            projindex = self.proj_index[ik, isp, ish, 0:dim]
        elif shells == 'all':
            if ir is None:
                raise ValueError("upfold: provide ir if treating all shells.")
            dim = self.shells[ish]['dim']
            # projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]
            projindex = self.proj_index_all[ik, isp, ish, ir, 0:dim]

        # gf_upfolded.from_L_G_R(
        #     projmat.conjugate().transpose(), gf_to_upfold, projmat)
        gf_upfolded.data[numpy.ix_(range(gf_to_upfold.data.shape[0]), projindex, projindex)] \
            += gf_to_upfold.data * fac

        if overwrite_gf_inp:
            return None
        else:
            return gf_upfolded

    ###############################################################
    # OVERRIDE FUNCTIONS
    # The density is computed by simple Matsubara summation
    # to avoid "High frequency error"
    # Replaced parts are indicated by "+++REPLACED"
    ###############################################################


    def total_density_matsubara(self, mu=None, iw_or_w="iw", with_Sigma=True, with_dc=True, broadening=None):
        r"""
        Calculates the total charge within the energy window for a given chemical potential.
        The chemical potential is either given by parameter `mu` or, if it is not specified,
        taken from `self.chemical_potential`.
        The total charge is calculated from the trace of the GF in the Bloch basis.
        By default, a full interacting GF is used. To use the non-interacting GF, set
        parameter `with_Sigma = False`.
        The number of bands within the energy windows generally depends on `k`. The trace is
        therefore calculated separately for each `k`-point.
        Since in general n_orbitals depends on k, the calculation is done in the following order:
        .. math:: n_{tot} = \sum_{k} n(k),
        with
        .. math:: n(k) = Tr G_{\nu\nu'}(k, i\omega_{n}).
        The calculation is done in the global coordinate system, if distinction is made between local/global.
        Parameters
        ----------
        mu : float, optional
             Input chemical potential. If not specified, `self.chemical_potential` is used instead.
        iw_or_w : string, optional
                  - `iw_or_w` = 'iw' for a imaginary-frequency self-energy
                  - `iw_or_w` = 'w' for a real-frequency self-energy
        with_Sigma : boolean, optional
             If `True` the full interacing GF is evaluated, otherwise the self-energy is not
             included and the charge would correspond to a non-interacting system.
        with_dc : boolean, optional
             Whether or not to subtract the double-counting term from the self-energy.
        broadening : float, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
                     Only relevant for real-frequency GF.
        Returns
        -------
        dens : float
               Total charge :math:`n_{tot}`.
        """

        if mu is None:
            mu = self.chemical_potential
        dens = 0.0
        ikarray = numpy.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):
            G_latt = self.lattice_gf(
                ik=ik, mu=mu, iw_or_w=iw_or_w, with_Sigma=with_Sigma, with_dc=with_dc, broadening=broadening)
            # dens += self.bz_weights[ik] * G_latt.total_density()
            # +++REPLACED
            dens += self.bz_weights[ik] * calc_total_density(G_latt)
        # collect data from mpi:
        dens = mpi.all_reduce(mpi.world, dens, lambda x, y: x + y)
        mpi.barrier()

        if abs(dens.imag) > 1e-20:
            mpi.report("Warning: Imaginary part in density will be ignored ({})".format(str(abs(dens.imag))))
        return dens.real


    def density_matrix_matsubara(self, method='using_gf', beta=40.0):
        """Calculate density matrices in one of two ways.
        Parameters
        ----------
        method : string, optional
                 - if 'using_gf': First get lattice gf (g_loc is not set up), then density matrix.
                                  It is useful for Hubbard I, and very quick.
                                  No assumption on the hopping structure is made (ie diagonal or not).
                 - if 'using_point_integration': Only works for diagonal hopping matrix (true in wien2k).
        beta : float, optional
               Inverse temperature.
        Returns
        -------
        dens_mat : list of dicts
                   Density matrix for each spin in each correlated shell.
        """
        dens_mat = [{} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            for sp in self.spin_block_names[self.corr_shells[icrsh]['SO']]:
                dens_mat[icrsh][sp] = numpy.zeros(
                    [self.corr_shells[icrsh]['dim'], self.corr_shells[icrsh]['dim']], numpy.complex_)

        ikarray = numpy.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):

            if method == "using_gf":

                G_latt_iw = self.lattice_gf(
                    ik=ik, mu=self.chemical_potential, iw_or_w="iw", beta=beta)
                # G_latt_iw *= self.bz_weights[ik]
                # dm = G_latt_iw.density()
                # MMat = [dm[sp] for sp in self.spin_block_names[self.SO]]
                # +++REPLACED
                dm = calc_density_matrix(G_latt_iw)
                MMat = [dm[sp] * self.bz_weights[ik] for sp in self.spin_block_names[self.SO]]

            elif method == "using_point_integration":

                ntoi = self.spin_names_to_ind[self.SO]
                spn = self.spin_block_names[self.SO]
                dims = {sp:self.n_orbitals[ik, ntoi[sp]] for sp in spn}
                MMat = [numpy.zeros([dims[sp], dims[sp]], numpy.complex_) for sp in spn]

                for isp, sp in enumerate(spn):
                    ind = ntoi[sp]
                    for inu in range(self.n_orbitals[ik, ind]):
                        # only works for diagonal hopping matrix (true in
                        # wien2k)
                        if (self.hopping[ik, ind, inu, inu] - self.h_field * (1 - 2 * isp)) < 0.0:
                            MMat[isp][inu, inu] = 1.0
                        else:
                            MMat[isp][inu, inu] = 0.0

            else:
                raise ValueError("density_matrix: the method '%s' is not supported." % method)

            for icrsh in range(self.n_corr_shells):
                for isp, sp in enumerate(self.spin_block_names[self.corr_shells[icrsh]['SO']]):
                    ind = self.spin_names_to_ind[
                        self.corr_shells[icrsh]['SO']][sp]
                    dim = self.corr_shells[icrsh]['dim']
                    n_orb = self.n_orbitals[ik, ind]
                    projmat = self.proj_mat[ik, ind, icrsh, 0:dim, 0:n_orb]
                    if method == "using_gf":
                        dens_mat[icrsh][sp] += numpy.dot(numpy.dot(projmat, MMat[isp]),
                                                         projmat.transpose().conjugate())
                    elif method == "using_point_integration":
                        dens_mat[icrsh][sp] += self.bz_weights[ik] * numpy.dot(numpy.dot(projmat, MMat[isp]),
                                                                               projmat.transpose().conjugate())

        # get data from nodes:
        for icrsh in range(self.n_corr_shells):
            for sp in dens_mat[icrsh]:
                dens_mat[icrsh][sp] = mpi.all_reduce(
                    mpi.world, dens_mat[icrsh][sp], lambda x, y: x + y)
        mpi.barrier()

        if self.symm_op != 0:
            dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for sp in dens_mat[icrsh]:
                    if self.rot_mat_time_inv[icrsh] == 1:
                        dens_mat[icrsh][sp] = dens_mat[icrsh][sp].conjugate()
                    dens_mat[icrsh][sp] = numpy.dot(numpy.dot(self.rot_mat[icrsh].conjugate().transpose(), dens_mat[icrsh][sp]),
                                                    self.rot_mat[icrsh])

        return dens_mat
