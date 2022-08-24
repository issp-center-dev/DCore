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
from collections import OrderedDict
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
        data[iw, i]  (i=(orb,spin))
    to
        g[spin].data[iw,orb1,orb2] (w/o SO) (spin='up', 'down')
        g['ud'].data[iw,i1,i2]     (w/  SO)

    """

    if g.n_blocks > 2:
        raise RuntimeError("n_blocks={} must be 1 or 2.".format(g.n_blocks))

    # (n_sp, n_inner) = (2, n_orb)    w/o SO
    #                   (1, 2*n_orb)  w/  SO
    n_sp = g.n_blocks
    n_inner = data.shape[1] // n_sp
    niw = data.shape[0]

    # array which data are assigned from
    # (orb, spin, iw) -> (spin, orb, iw) -> (spin, orb, iw)  w/o SO
    #                                       (1, inner, iw)   w/  SO
    # data_from = data.reshape((-1, 2, niw)).transpose((1, 0, 2)).reshape(n_sp, n_inner, niw)

    # (iw, orb, spin) -> (iw, spin, orb) -> (iw, spin, orb)  w/o SO
    #                                       (iw, 1, inner)   w/  SO
    data_from = data.reshape((niw, -1, 2)).transpose((0, 2, 1)).reshape(niw, n_sp, n_inner)

    # assign data_from to g
    for sp, name in enumerate(names):
        # array which data are assigned to
        data_to = g[name].data
        assert data_to.shape == (2*niw, n_inner, n_inner)
        for i in range(n_inner):
            # positive frequency
            data_to[niw:, i, i] = data_from[:, sp, i]
            # negative frequency
            data_to[:niw, i, i] = numpy.conj(data_from[::-1, sp, i])


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

def convert_Umatrix(U, Uprime, J, norb):
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

    return Uout.reshape((2*norb, 2*norb))

    # Uout = Uout.reshape((2*norb, 2*norb))
    # with open('./Umatrix', 'w') as f:
    #     for i in range(2*norb):
    #         for j in range(2*norb):
    #             print('{:.15e} '.format(Uout[i, j].real), file=f, end="")
    #         print("", file=f)


def eval_tail(g_iw, beta, n_ave=10):
    # array_iw[iw, flavor]
    assert g_iw.ndim == 2

    n_iw, n_flavors = g_iw.shape
    iw = (2*numpy.arange(n_iw) + 1) * numpy.pi / beta

    tails = numpy.mean(g_iw[-n_ave:, :] * iw[-n_ave:, None], axis=0).real
    assert tails.shape == (n_flavors,)

    return tails


class JOCTHYBSEGSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.
        """

        super().__init__(beta, gf_struct, u_mat, n_iw)

        # self.n_tau = max(10001, 5 * n_iw)

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

        # (0) Rotate H0 and Delta_tau if rot is given
        g0_iw_rotated = self._G0_iw.copy()
        if rot is None:
            u_mat_rotated = self.u_mat
        else:
            assert isinstance(rot, dict)
            u_mat_rotated = rotate_basis(rot, self.use_spin_orbit, self.u_mat, [g0_iw_rotated,], direction='forward')

        # (1) Set configuration for the impurity solver
        # input:
        #   self.beta
        #   self.set_G0_iw
        #   self.u_mat
        #
        # Additionally, the following variables may used:
        #   self.n_orb
        #   self.n_flavors
        #   self.gf_struct

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
        # H0 = (H0[:, index])[index, :]
        H0 = H0[numpy.ix_(index, index)]

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

        with open('./ef.in', 'w') as f:
            for i in range(2*self.n_orb):
                print('{:.15e}'.format(H0[i, i].real), file=f)


        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)

        # H0 is extracted from the tail of the Green's function.
        self._Delta_iw = delta(g0_iw_rotated)
        # Delta_tau = make_block_gf(GfImTime, self.gf_struct, self.beta, self.n_tau)
        # for name in self.block_names:
        #     Delta_tau[name] << Fourier(self._Delta_iw[name])
        # Delta_tau_data = to_numpy_array(Delta_tau, self.block_names)
        delta_iw = to_numpy_array(self._Delta_iw, self.block_names)
        print(delta_iw.shape)
        assert delta_iw.shape == (self.n_iw * 2, self.n_orb * 2, self.n_orb * 2)
        delta_iw = delta_iw[self.n_iw:, :, :]  # only positive frequency
        print(delta_iw.shape)
        if self.n_iw % 2 != 0:
            sys.exit(f"Invalid value n_iw={self.n_iw}: Only even number is allowed for n_iw in JO/cthyb-seg solver")

        # TODO: check delta_iw
        #    Delta_{ab}(iw) must be diagonal

        print(delta_iw.shape)
        delta_iw_diagonal = numpy.einsum("wii->wi", delta_iw)
        print(delta_iw_diagonal.shape)

        with open('./delta.in', 'w') as f:
            for iw in range(self.n_iw):
                for f1 in range(self.n_flavors):
                    val = delta_iw_diagonal[iw, f1]
                    print(' {:.15e} {:.15e}'.format(val.real, val.imag), file=f, end="")
                print("", file=f)

        # tail of Delta(iw)
        vsq = eval_tail(delta_iw_diagonal, self.beta)
        with open('./vsq.in', 'w') as f:
            for i in range(self.n_flavors):
                print('{:.15e}'.format(vsq[i]), file=f)


        # (1c) Set U_{ijkl} for the solver
        # Set up input parameters and files for ALPS/CTHYB-SEG

        U, Uprime, J = dcore2alpscore(u_mat_rotated)
        # write_Umatrix(U, Uprime, J, self.n_orb)
        Udd = convert_Umatrix(U, Uprime, J, self.n_orb)
        assert Udd.shape == (self.n_flavors, self.n_flavors)

        with open('./u.in', 'w') as f:
            for i in range(self.n_flavors):
                for j in range(self.n_flavors):
                    print(' {:.15e}'.format(Udd[i, j].real), file=f, end="")
                print("", file=f)


        # (1d) Set parameters for the solver

        params_solver = OrderedDict()
        params_solver['model'] = {
            'n_s' : self.n_orb*2,
            'file_Delta' : 'delta.in',
            'file_Vsq' : 'vsq.in',
            'file_U' : 'u.in',
            'file_ef' : 'ef.in',
            'beta' : self.beta,
        }
        params_solver['control'] = {
            'n_tau' : self.n_iw * 2,
        }
        params_solver['MC'] = {}

        # Set parameters specified by the user
        for prefix in ['control', 'MC']:
            for k, v in params_kw.items():
                # Convert, e.g., MC.n_msr -> n_msr
                if k.startswith(prefix + '.'):
                    key = k[len(prefix)+1:]
                    if key in params_solver[prefix]:
                        sys.exit(f"ERROR: Cannot override parameter '{key}'")
                    params_solver[prefix][key] = v

        with open('./input.ini', 'w') as f:
            for prefix, params in params_solver.items():
                print(f"\n[{prefix}]", file=f)
                for k, v in params.items():
                    print(f"{k} = {v}", file=f)

        if _read('dry_run'):
            return

        # (2) Run a working horse
        exec_path = expand_path(_read('exec_path'))
        with open('./output', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, './input.ini'], output_f)

        with open('./output', 'r') as output_f:
            for line in output_f:
                print(line, end='')

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

        data = numpy.loadtxt("self_w.dat")
        print(data.shape)
        sigma_data = data[:, 1::2] + 1j * data[:, 2::2]
        print(sigma_data.shape)
        assert sigma_data.shape == (self.n_iw, self.n_flavors)
        assign_from_numpy_array(self._Sigma_iw, sigma_data, self.block_names)

        data = numpy.loadtxt("Gf_w.dat")
        print(data.shape)
        gf_data = data[:, 1::2] + 1j * data[:, 2::2]
        print(gf_data.shape)
        assert gf_data.shape == (self.n_iw, self.n_flavors)
        assign_from_numpy_array(self._Gimp_iw, gf_data, self.block_names)

        # Rotate Sigma and Gimp back to the original basis
        if rot is not None:
            rotate_basis(rot, self.use_spin_orbit, None, [self._Sigma_iw, self._Gimp_iw], direction='backward')
            # rotate_basis(rot, self.use_spin_orbit, None, self._Gimp_iw, direction='backward')

        #   self.quant_to_save['nn_equal_time']
        # nn_equal_time = self._get_results("nn", 1, orbital_symmetrize=True, stop_if_data_not_exist=False)
        # [(s1,o1), (s2,o2), 0]
        # self.quant_to_save['nn_equal_time'] = nn_equal_time[:, :, 0]  # copy

    def calc_Xloc_ph(self, rot, mpirun_command, num_wf, num_wb, params_kw):
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
        raise NotImplementedError

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
        return "JO/cthyb-seg"
