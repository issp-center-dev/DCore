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
import os
import shutil
import sys
from itertools import product
from dcore._dispatcher import *
from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0, umat2dd, get_block_size, expand_path
from .base import SolverBase, rotate_basis


def to_numpy_array(g, names):
    """
    Convert BlockGf object to numpy.
    Rearrange spins and orbitals so that up and down spins appear alternatively.
    If there is a single block, we assume that spin and down spins appear alternatively.
    If there are two blocks, we assume that they are spin1 and spin2 sectors.
    """

    if g.n_blocks > 2:
        raise RuntimeError("n_blocks={} must be 1 or 2.".format(g.n_blocks))

    n_spin_orbital = numpy.sum([get_block_size(block) for _, block in g])

    # FIXME: Bit ugly
    n_data = g[names[0]].data.shape[0]

    data = numpy.zeros((n_data, n_spin_orbital, n_spin_orbital), dtype=complex)
    offset = 0
    for name in names:
        block = g[name]
        block_dim = get_block_size(block)
        data[:, offset:offset + block_dim, offset:offset + block_dim] = block.data
        offset += block_dim

    # from (spin, orb) : (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
    # to (orb, spin)   : (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
    norb = n_spin_orbital//2
    index = numpy.zeros(n_spin_orbital, dtype=int)
    index[0::2] = numpy.arange(norb)
    index[1::2] = numpy.arange(norb) + norb
    # Swap cols and rows
    return (data[:, :, index])[:, index, :]


# def assign_from_numpy_array(g, data, names):
#     """
#     Does inversion of to_numpy_array
#     data[spin,orb,iw]
#     g[spin].data[iw,orb1,orb2]
#     """
#     # print(g.n_blocks)
#     if g.n_blocks != 2:
#         raise RuntimeError("n_blocks={} must be 1 or 2.".format(g.n_blocks))
#
#     norb = data.shape[1]
#     niw = data.shape[2]
#     # print(data.shape)
#     # print("norb:", norb)
#
#     # check number of Matsubara frequency
#     assert data.shape[2]*2 == g[names[0]].data.shape[0]
#     # print(g[names[0]].data.shape)
#
#     for spin in range(2):
#         for orb in range(norb):
#             # print(orb, spin, names[spin])
#             # positive frequency
#             g[names[spin]].data[niw:, orb, orb] = data[spin][orb][:]
#             # negative frequency
#             g[names[spin]].data[:niw, orb, orb] = numpy.conj(data[spin][orb][::-1])


def assign_from_numpy_array(g, data, names):
    """
    Does inversion of to_numpy_array

    assign from
        data[i, iw]  (i=(orb,spin))
    to
        g[spin].data[iw,orb1,orb2] (w/o SO) (spin='up', 'down')
        g['ud'].data[iw,i1,i2]     (w/  SO)

    """

    if g.n_blocks > 2:
        raise RuntimeError("n_blocks={} must be 1 or 2.".format(g.n_blocks))

    # (n_sp, n_inner) = (2, n_orb)    w/o SO
    #                   (1, 2*n_orb)  w/  SO
    n_sp = g.n_blocks
    n_inner = data.shape[0] // n_sp
    niw = data.shape[1]

    # array which data are assigned from
    # (orb, spin, iw) -> (spin, orb, iw) -> (spin, orb, iw)  w/o SO
    #                                       (1, inner, iw)   w/  SO
    data_from = data.reshape((-1, 2, niw)).transpose((1, 0, 2)).reshape(n_sp, n_inner, niw)

    # assign data_from to g
    for sp, name in enumerate(names):
        # array which data are assigned to
        data_to = g[name].data
        assert data_to.shape == (2*niw, n_inner, n_inner)
        for i in range(n_inner):
            # positive frequency
            data_to[niw:, i, i] = data_from[sp, i, :]
            # negative frequency
            data_to[:niw, i, i] = numpy.conj(data_from[sp, i, ::-1])


def dcore2alpscore(dcore_U):

    dcore_U_len = len(dcore_U)
    alps_U = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_Uprime = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)
    alps_J = numpy.zeros((dcore_U_len, dcore_U_len), dtype=float)

    # m_range = range(size)
    for i, j in product(range(dcore_U_len), range(dcore_U_len)):
        alps_U[i, j] = dcore_U[i, j, i, j].real - dcore_U[i, j, j, i].real
        alps_Uprime[i, j] = dcore_U[i, j, i, j].real
        alps_J[i, j] = dcore_U[i, j, j, i].real
    return alps_U, alps_Uprime, alps_J

def write_Umatrix(U, Uprime, J, norb):
    Uout = numpy.zeros((norb, 2, norb, 2))

    # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
    # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
    def func(u):
        uout = u.reshape((2, norb, 2, norb)).transpose(1, 0, 3, 2)
        return uout

    U_four = func(U)
    Uprime_four = func(Uprime)
    J_four = func(J)

    for a1, a2 in product(range(norb), repeat=2):
        for s1, s2 in product(range(2), repeat=2):  # spin-1/2
            if a1 == a2:
                Uout[a1, s1, a2, s2] = U_four[a1, s1, a2, s2]
            else:
                Uout[a1, s1, a2, s2] = Uprime_four[a1, s1, a2, s2] - J_four[a1, s1, a2, s2]

    Uout = Uout.reshape((2*norb, 2*norb))
    with open('./Umatrix', 'w') as f:
        for i in range(2*norb):
            for j in range(2*norb):
                print('{:.15e} '.format(Uout[i, j].real), file=f, end="")
            print("", file=f)



class ALPSCTHYBSEGSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.
        """

        super(ALPSCTHYBSEGSolver, self).__init__(beta, gf_struct, u_mat, n_iw)

        self.n_tau = max(10001, 5 * n_iw)

    def _get_occupation(self):
        """
        Read the spin-orbital-dependent occupation number from HDF5 file

        Returns
        -------
        numpy.ndarray of size (2*self.n_orb)

        """

        array = numpy.zeros(2*self.n_orb, dtype=float)
        with HDFArchive('sim.h5', 'r') as f:
            results = f["simulation"]["results"]
            for i1 in range(2*self.n_orb):
                group = "density_%d" % i1
                if group in results:
                    array[i1] = results[group]["mean"]["value"]

        # [(o1,s1)] -> [o1, s1] -> [s1, o1] -> [(s1, o1)]
        array = array.reshape((self.n_orb, 2))\
                     .transpose((1, 0))\
                     .reshape((2*self.n_orb))
        return array

    def _get_results(self, group_prefix, n_data, orbital_symmetrize, dtype=float, stop_if_data_not_exist=True):
        """
        Read results with two spin-orbital indices from HDF5 file

        Returns
        -------
        numpy.ndarray of size (2*self.n_orb, 2*self.n_orb, n_data)

        """

        data_shape = (2*self.n_orb, 2*self.n_orb, n_data)

        array = numpy.zeros(data_shape, dtype=dtype)
        with HDFArchive('sim.h5', 'r') as f:
            results = f["simulation"]["results"]
            for i1, i2 in product(range(2*self.n_orb), repeat=2):
                group = "%s_%d_%d" % (group_prefix, i1, i2)
                if group in results:
                    array[i1, i2, :] = results[group]["mean"]["value"]
                    if orbital_symmetrize:  # Only i1>i2 is computed in CTQMC.
                        array[i2, i1, :] = array[i1, i2, :]
                elif stop_if_data_not_exist:
                    raise Exception("data does not exist in sim.h5/simulation/results/{}. alps_cthyb might be old.".format(group))


        # [(o1,s1), (o2,s2)] -> [o1, s1, o2, s2] -> [s1, o1, s2, o2] -> [(s1,o1), (s2,o2)]
        array = array.reshape((self.n_orb, 2, self.n_orb, 2, -1))\
                     .transpose((1, 0, 3, 2, 4))\
                     .reshape((2*self.n_orb, 2*self.n_orb, -1))
        return array

    def solve(self, rot, mpirun_command, params_kw):
        """
        In addition to the parameters described in the docstring of SolverBase,
        one can pass solver-dependent parameters using params_kw. For example,
          exec_path : str, path to an executable, mandatory
          dry_run   : bool, actual computation is not performed if dry_run is True, optional
        """
        internal_params = {
            'exec_path'           : '',
            'random_seed_offset'  : 0,
            'dry_run'             : False,
            'neglect_offdiagonal' : True,
        }

        def _read(key):
            if key in params_kw:
                return params_kw[key]
            else:
                return internal_params[key]
        print (params_kw)

        umat_check = umat2dd(self.u_mat)
        assert numpy.allclose(umat_check, self.u_mat), "Please set density_density = True when you run ALPS/cthyb-seg!"

        # (1) Set configuration for the impurity solver
        # input:
        #   self.beta
        #   self.set_G0_iw
        #   self.u_mat
        #
        # Additionally, the following variables may used:
        #   self.n_orb
        #   self.n_flavor
        #   self.gf_struct

        # (0) Rotate H0 and Delta_tau if rot is given
        g0_iw_rotated = self._G0_iw.copy()
        if rot is None:
            u_mat_rotated = self.u_mat
        else:
            assert isinstance(rot, dict)
            u_mat_rotated = rotate_basis(rot, self.use_spin_orbit, self.u_mat, [g0_iw_rotated,], direction='forward')

        # (1a) If H0 is necessary:
        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin1, ..., spin2, spin2, ...
        H0 = extract_H0(g0_iw_rotated, self.block_names)

        # from (up,orb1), (up,orb2), ..., (down,orb1), (down,orb2), ...
        # to (up,orb1), (down,orb1), (up,orb2), (down,orb2), ...
        index = numpy.zeros((2*self.n_orb), dtype=int)
        index[0::2] = numpy.arange(self.n_orb)
        index[1::2] = numpy.arange(self.n_orb) + self.n_orb
        # Swap cols and rows
        H0 = (H0[:, index])[index, :]

        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)
        # H0 is extracted from the tail of the Green's function.
        self._Delta_iw = delta(g0_iw_rotated)
        Delta_tau = make_block_gf(GfImTime, self.gf_struct, self.beta, self.n_tau)
        for name in self.block_names:
            Delta_tau[name] << Fourier(self._Delta_iw[name])
        Delta_tau_data = to_numpy_array(Delta_tau, self.block_names)

        # check if H0 is diagonal
        H0_offdiag = H0.copy()
        for i in range(H0_offdiag.shape[0]):
            H0_offdiag[i, i] = 0
        if numpy.linalg.norm(H0_offdiag) > 1e-6:
            print("\nWARNING: The local Hamiltonian is not diagonal", file=sys.stderr)
            print("H0_loc =\n{}".format(H0), file=sys.stderr)
            if _read('neglect_offdiagonal'):
                print("--> continue. To stop calculation, set neglect_offdiagonal{bool}=False", file=sys.stderr)
            else:
                print("--> exit. To neglect this warning, set neglect_offdiagonal{bool}=True", file=sys.stderr)
                sys.exit(1)

        # TODO: check Delta_tau_data
        #    Delta_{ab}(tau) should be diagonal, real, negative

        # (1c) Set U_{ijkl} for the solver
        # Set up input parameters and files for ALPS/CTHYB-SEG

        p_run = {
            'SEED'                            : params_kw['random_seed_offset'],
            'FLAVORS'                         : self.n_orb*2,
            'BETA'                            : self.beta,
            'N'                               : self.n_tau - 1,
            'NMATSUBARA'                      : self.n_iw,
            'U_MATRIX'                        : 'Umatrix',
            'MU_VECTOR'                       : 'MUvector',
            'cthyb.DELTA'                     : 'delta',
        }

        if os.path.exists('./input.out.h5'):
            shutil.move('./input.out.h5', './input_prev.out.h5')
        # Set parameters specified by the user
        for k, v in params_kw.items():
            if k in internal_params:
                continue
            if k in p_run:
                raise RuntimeError("Cannot override input parameter for ALPS/CTHYB-SEG: " + k)
            p_run[k] = v

        with open('./input.ini', 'w') as f:
            for k, v in p_run.items():
                print(k, " = ", v, file=f)

        with open('./delta', 'w') as f:
            for itau in range(self.n_tau):
                print('{}'.format(itau), file=f, end="")
                for f1 in range(self.n_flavors):
                    if Delta_tau_data[itau, f1, f1].real >0:
                        Delta_tau_data[itau, f1, f1] = 0
                    print(' {:.15e}'.format(Delta_tau_data[itau, f1, f1].real), file=f, end="")
                print("", file=f)

        U, Uprime, J = dcore2alpscore(u_mat_rotated)
        write_Umatrix(U, Uprime, J, self.n_orb)

        with open('./MUvector', 'w') as f:
            # for orb in range(self.n_orb):
            #     for spin in range(2):
            #         print('{:.15e} '.format(-H0[2*orb+spin][2*orb+spin].real), file=f, end="")
            for i in range(2*self.n_orb):
                print('{:.15e} '.format(-H0[i, i].real), file=f, end="")
            print("", file=f)

        if _read('dry_run'):
            return

        # Invoke subprocess
        exec_path = expand_path(_read('exec_path'))

        # (2) Run a working horse
        with open('./output', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, './input.ini'], output_f)

        with open('./output', 'r') as output_f:
            for line in output_f:
                print(line, end='')

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

        def set_blockgf_from_h5(sigma, group):
            # swdata = numpy.zeros((2, self.n_orb, self.n_iw), dtype=complex)
            swdata = numpy.zeros((2*self.n_orb, self.n_iw), dtype=complex)
            with HDFArchive('sim.h5', 'r') as f:
                # for orb in range(self.n_orb):
                #     for spin in range(2):
                #         swdata_array = f[group][str(orb*2+spin)]["mean"]["value"]
                #         assert swdata_array.dtype == numpy.complex
                #         assert swdata_array.shape == (self.n_iw,)
                #         swdata[spin, orb, :] = swdata_array
                for i in range(2*self.n_orb):
                    swdata_array = f[group][str(i)]["mean"]["value"]
                    assert swdata_array.dtype == numpy.complex
                    assert swdata_array.shape == (self.n_iw,)
                    swdata[i, :] = swdata_array
            assign_from_numpy_array(sigma, swdata, self.block_names)

        set_blockgf_from_h5(self._Sigma_iw, "S_omega")
        set_blockgf_from_h5(self._Gimp_iw, "G_omega")

        # Rotate Sigma and Gimp back to the original basis
        if rot is not None:
            rotate_basis(rot, self.use_spin_orbit, None, [self._Sigma_iw, self._Gimp_iw], direction='backward')
            # rotate_basis(rot, self.use_spin_orbit, None, self._Gimp_iw, direction='backward')

        #   self.quant_to_save['nn_equal_time']
        nn_equal_time = self._get_results("nn", 1, orbital_symmetrize=True, stop_if_data_not_exist=False)
        # [(s1,o1), (s2,o2), 0]
        self.quant_to_save['nn_equal_time'] = nn_equal_time[:, :, 0]  # copy

    def calc_G2loc_ph(self, rot, mpirun_command, num_wf, num_wb, params_kw):
        """
        compute local G2 in p-h channel
            X_loc = < c_{i1}^+ ; c_{i2} ; c_{i4}^+ ; c_{i3} >

        Parameters
        ----------
        rot
        mpirun_command
        num_wf
        num_wb
        params_kw

        Returns
        -------
        G2_loc : dict
            key = (i1, i2, i3, i4)
            val = numpy.ndarray(n_w2b, 2*n_w2f, 2*n_w2f)

        chi_loc : dict (None if not computed)
            key = (i1, i2, i3, i4)
            val = numpy.ndarray(n_w2b)
        """
        if rot is not None:
            # TODO
            raise NotImplementedError

        use_chi_loc = False

        params_kw['cthyb.MEASURE_g2w'] = 1
        params_kw['cthyb.N_w2'] = num_wf
        params_kw['cthyb.N_W'] = num_wb
        if use_chi_loc:
            params_kw['cthyb.MEASURE_nnw'] = 1

        self.solve(rot, mpirun_command, params_kw)

        # Save G2(wb, wf, wf')
        # [(s1,o1), (s2,o2), (wb,wf,wf')]
        g2_re = self._get_results("g2w_re", 4*num_wf*num_wf*num_wb, orbital_symmetrize=False)
        g2_im = self._get_results("g2w_im", 4*num_wf*num_wf*num_wb, orbital_symmetrize=False)
        g2_loc = (g2_re + g2_im * 1.0J) / self.beta
        g2_loc = g2_loc.reshape((2*self.n_orb, 2*self.n_orb) + (num_wb, 2*num_wf, 2*num_wf))
        # assign to dict
        g2_dict = {}
        for i1, i2 in product(range(2*self.n_orb), repeat=2):
            g2_dict[(i1, i1, i2, i2)] = g2_loc[i1, i2]

        # return g2_loc for arbitrary wb including wb<0
        def get_g2(_i, _j, _wb, _wf1, _wf2):
            try:
                if _wb >= 0:
                    return g2_loc[_i, _j, _wb, _wf1, _wf2]
                else:
                    # G2_iijj(wb, wf, wf') = G2_jjii(-wb, -wf', -wf)^*
                    return numpy.conj(g2_loc[_j, _i, -_wb, -(1+_wf2), -(1+_wf1)])
            except IndexError:
                return 0

        # Convert G2_iijj -> G2_ijij
        g2_loc_tr = numpy.zeros(g2_loc.shape, dtype=complex)
        for i1, i2 in product(range(2*self.n_orb), repeat=2):
            for wb in range(num_wb):
                for wf1, wf2 in product(range(2 * num_wf), repeat=2):
                    # G2_ijij(wb, wf, wf') = -G2_iijj(wf-wf', wf'+wb, wf')^*
                    g2_loc_tr[i1, i2, wb, wf1, wf2] = -get_g2(i1, i2, wf1-wf2, wf2+wb, wf2)
        # assign to dict
        for i1, i2 in product(range(2*self.n_orb), repeat=2):
            # exclude i1=i2, which was already assigned by g2_loc
            if i1 != i2:
                g2_dict[(i1, i2, i1, i2)] = g2_loc_tr[i1, i2]

        # Occupation number
        # [(s1,o1)]
        occup = self._get_occupation()

        # Save chi(wb)
        # [(s1,o1), (s2,o2), wb]
        chi_dict = None
        if use_chi_loc:
            chi_re = self._get_results("nnw_re", num_wb, orbital_symmetrize=True)
            chi_im = self._get_results("nnw_im", num_wb, orbital_symmetrize=True)
            chi_loc = chi_re + chi_im * 1.0J
            # subtract <n><n>
            chi_loc[:, :, 0] -= occup[:, None] * occup[None, :] * self.beta
            # assign to dict
            chi_dict = {}
            for i1, i2 in product(range(2*self.n_orb), repeat=2):
                chi_dict[(i1, i1, i2, i2)] = chi_loc[i1, i2]

        return g2_dict, chi_dict

    def calc_G2loc_ph_sparse(self, rot, mpirun_command, freqs_ph, num_wb, params_kw):
        raise Exception("This solver does not support the sparse sampling.")

    def name(self):
        return "ALPS/cthyb-seg"
