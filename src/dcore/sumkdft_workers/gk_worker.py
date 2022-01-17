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

from .._dispatcher import mpi

from ..mpi import split_idx, gatherv
from ..sumkdft_opt import SumkDFT_opt
from ..tools import gf_block_names
from .worker_base import SumkDFTWorkerBase, setup_sk

def _from_L_G_R(L, G, R):
    return numpy.einsum('ij,wjk,kl->wil', L, G, R)

class SumkDFTGk(SumkDFT_opt):
    """
    Extends the SumkDFT class for computing intra/inter shells components of G(k)
    """
    def __init__(self, hdf_file, **kwargs):
        """
        Initialisation of the class. Parameters are exactly as for SumKDFT.
        """
        super().__init__(hdf_file, **kwargs)


    def downfold_offdiagonal(self, ik, ish1, ish2, bname, gf_to_downfold, shells='corr', ir=None):
        r"""
        [ORIGINAL: function 'downfold' in sumk_dft.py]
        Downfolds a block of the Green's function for a given shell and k-point using the corresponding projector matrices.
        Returns the result a ndarray object.

        Parameters
        ----------
        ik : integer
             k-point index for which the downfolding is to be done.
        ish : integer
              Shell index of GF to be downfolded.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        bname : string
                Block name of the target block of the lattice Green's function.
        gf_to_downfold : Gf
                       Block of the Green's function that is to be downfolded.
        shells : string, optional

                 - if shells='corr': orthonormalized projectors for correlated shells are used for the downfolding.
                 - if shells='all': non-normalized projectors for all included shells are used for the downfolding.

        ir : integer, optional
             Index of equivalent site in the non-correlated shell 'ish', only used if shells='all'.

        Returns
        -------
        gf_downfolded : numpy.ndarray (npoints, nso, nso)
                      Downfolded block of the lattice Green's function.
        """

        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim1 = self.corr_shells[ish1]['dim']
            dim2 = self.corr_shells[ish2]['dim']
            projmat1 = self.proj_mat[ik, isp, ish1, 0:dim1, 0:n_orb]
            projmat2 = self.proj_mat[ik, isp, ish2, 0:dim2, 0:n_orb]
        elif shells == 'all':
            if ir is None:
                raise ValueError("downfold: provide ir if treating all shells.")
            dim1 = self.shells[ish1]['dim']
            dim2 = self.shells[ish2]['dim']
            projmat1 = self.proj_mat_all[ik, isp, ish1, ir, 0:dim1, 0:n_orb]
            projmat2 = self.proj_mat_all[ik, isp, ish2, ir, 0:dim2, 0:n_orb]
        
        return _from_L_G_R(projmat1, gf_to_downfold.data, projmat2.conjugate().transpose())


    def extract_Gk_crsh_sparse(self, smpl_freqs, mu=None, iw_or_w='iw', with_Sigma=True, with_dc=True, broadening=None):
        r"""
        [ORIGINAL: function 'extract_G_loc' in sumk_dft.py]
        Extracts the k-resolved downfolded Green function.

        Parameters
        ----------
        mu : real, optional
             Input chemical potential. If not provided the value of self.chemical_potential is used as mu.
        smpl_freqs : 1D int array
             Indices of sampling fermionic frequencies (e.g., -1, 0, 1, ...)
        with_Sigma : boolean, optional
                     If True then the local GF is calculated with the self-energy self.Sigma_imp.
        with_dc : boolean, optional
                  If True then the double-counting correction is subtracted from the self-energy in calculating the GF.
        broadening : float, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
                     Only relevant for real-frequency GF.

        Returns
        -------
        G_k : numpy.ndarray of shape (nk_local, smpl_freqs.size, num_so, num_so)
            Gk rotated into the corresponding local frames, distributed over MPI processes
        """
        assert broadening is None
        if iw_or_w != "iw":
            raise Exception("Only iw_or_w=='iw' implemented")
        if mu is None:
            mu = self.chemical_potential

        # icrsh=0
        mesh = self.Sigma_imp_iw[0].mesh
        beta = mesh.beta
        use_spin_orbit = self.Sigma_imp_iw[0].n_blocks == 1
        block_names = gf_block_names(use_spin_orbit)
        niw = len(mesh)//2 # Number of non-negative frequencies
        if mesh.positive_only():
            raise RuntimeError("mesh must not be positive_only!")

        # Smpl frequencies
        if numpy.abs(smpl_freqs).max() >= niw:
            raise RuntimeError("Some of sampling frequencies are not on the frequency mesh. "
                       "Data on these sampling points will be replaced by 1/iw."
                       "This may cause systematic errors if the frequency window size is small."
                       )
        smpl_freqs_idx = smpl_freqs + niw

        # Number of spin orbitals in a unit cell
        # dim_corr_sh = 2 * num_orb (use_spin_orbit=True)
        #             =     num_orb (use_spin_orbit=False)
        dim_corr_sh = numpy.array([sh['dim'] for sh in self.corr_shells])
        num_so_corr_sh = dim_corr_sh if use_spin_orbit else 2*dim_corr_sh
        num_orb_corr_sh = num_so_corr_sh//2
        num_so = numpy.sum(num_so_corr_sh)

        # Distribute k points over MPI processes
        local_ksizes, k_offsets = split_idx(self.n_k, mpi.size)

        _offsets = numpy.hstack((0, numpy.cumsum(num_orb_corr_sh)))
        gk_local = numpy.zeros((local_ksizes[mpi.rank], smpl_freqs.size, 2, num_so//2, 2, num_so//2), dtype=numpy.complex128)
        for ik_local in range(local_ksizes[mpi.rank]):
            ik = ik_local + k_offsets[mpi.rank]
            G_latt = self.lattice_gf(ik=ik,
                mu=mu, iw_or_w=iw_or_w, with_Sigma=with_Sigma, with_dc=with_dc, beta=beta)
            for icrsh1 in range(self.n_corr_shells):
                for icrsh2 in range(self.n_corr_shells):
                    norb_sh1 = num_orb_corr_sh[icrsh1]
                    norb_sh2 = num_orb_corr_sh[icrsh2]
                    if use_spin_orbit:
                        gf = self.downfold_offdiagonal(ik, icrsh1, icrsh2, 'ud', G_latt['ud'])
                        if self.use_rotations:
                            # Rotate to the local coordinate system:
                            gf = _from_L_G_R(self.rot_mat[icrsh1].conjugate().T, gf, self.rot_mat[icrsh2])
                        gf = gf.reshape((gf.shape[0], 2, norb_sh1, 2, norb_sh2))
                        gk_local[ik_local, :,
                                :, _offsets[icrsh1]:_offsets[icrsh1+1], # spin, orb
                                :, _offsets[icrsh2]:_offsets[icrsh2+1] # spin, orb
                            ] = gf[smpl_freqs_idx, :, :, :, :]
                    else:
                        for ib, bname in enumerate(block_names): # 'up', 'down'
                            gf = self.downfold_offdiagonal(ik, icrsh1, icrsh2, bname, G_latt[bname])
                            if self.use_rotations:
                                # Rotate to the local coordinate system:
                                gf = _from_L_G_R(self.rot_mat[icrsh1].conjugate().T, gf, self.rot_mat[icrsh2])
                            gf = gf.reshape((gf.shape[0], norb_sh1, norb_sh2))
                            gk_local[ik_local, :,
                                    ib, _offsets[icrsh1]:_offsets[icrsh1+1], # spin, orb
                                    ib, _offsets[icrsh2]:_offsets[icrsh2+1]  # spin, orb
                                ] = gf[smpl_freqs_idx, :, :]
        return gk_local


class SumkDFTWorkerGk(SumkDFTWorkerBase):
    """For computing Gk"""
    def __init__(self, model_hdf5_file, input_file, output_file) -> None:
        super().__init__(model_hdf5_file, input_file, output_file)

    def run(self):
        with_dc = self.params['with_dc']
        smpl_freqs = self.params['smpl_freqs']
        nfreqs = len(smpl_freqs)

        sk = SumkDFTGk(self.model_hdf5_file)
        setup_sk(sk, 'iwn', self.params)

        # k-dependent Green's function
        gk_local = sk.extract_Gk_crsh_sparse(smpl_freqs, with_dc=with_dc)

        # gk_local: (local_ksizes[mpi.rank], smpl_freqs.size, 2, num_orb, 2, num_orb)
        num_orb = gk_local.shape[3]
        num_so = 2 * num_orb
        gk = gatherv(gk_local.ravel(), root=0)

        if mpi.is_master_node():
            gk = gk.reshape((sk.n_k, nfreqs, 2, num_so//2, 2, num_so//2))
            self.save_result({'gk': gk})