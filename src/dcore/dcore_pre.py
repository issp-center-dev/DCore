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
import sys
import os
import numpy
import ast
import h5py
from itertools import product
from dcore._dispatcher import HDFArchive, U_J_to_radial_integrals, U_matrix, cubic_names
from dcore.program_options import create_parser

from dcore.converters.wannier90 import Wannier90Converter

from dcore.tools import *
from dcore.sumkdft_compat import SumkDFTCompat

from dcore.lattice_models import create_lattice_model
from dcore.lattice_models.tools import print_local_fields
from dcore.program_options import parse_parameters


def __print_paramter(p, param_name):
    print(param_name + " = " + str(p[param_name]))


def _check_parameters(p, required, unused):
    for key in required:
        if p[key] == "None":
            print(f"Error ! Parameter '{key}' is not specified.")
            sys.exit(-1)
    for key in unused:
        if p[key] != "None":
            print(f"Error ! Parameter '{key}' is specified but is not used.")
            sys.exit(-1)


def _parse_interaction_parameters(input_str, name, nsh, n_inner=None):
    # parse interaction parameters
    # return list of list
    try:
        list_of_list = ast.literal_eval(input_str)
        if not isinstance(list_of_list, (list, tuple)):
            raise Exception(f"type({name}) must be list or tuple but is {type(list_of_list)}")
        if len(list_of_list) != nsh:
            raise Exception(f"len({name}) must be {nsh} but is {len(list_of_list)}")
        for uvalues in list_of_list:
            if not isinstance(uvalues, (list, tuple)):
                raise Exception(f"type({uvalues!r}) must be list or tuple but is {type(uvalues)}")
            if n_inner is not None:
                if len(uvalues) != n_inner:
                    raise Exception(f"len({uvalues!r}) must be {n_inner} but is {len(uvalues)}")
    except Exception as e:
        print(f"\nERROR in parsing {name} = {input_str!r}", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(-1)
    print(f" {name} = {list_of_list!r}")
    return list_of_list


def _generate_umat_kanamori(p):
    _check_parameters(p['model'], required=['kanamori'], unused=['slater_f', 'slater_uj'])

    # parse kanamori parameters
    nsh = p['model']['n_inequiv_shells']
    kanamori_sh = _parse_interaction_parameters(p["model"]["kanamori"], name='kanamori', nsh=nsh, n_inner=3)

    # print summary
    print("\n Kanamori interactions")
    for ish, (_u, _up, _j) in enumerate(kanamori_sh):
        print(f"  ish = {ish}")
        print(f"    | U = {_u!r}")
        print(f"    | U' = {_up!r}")
        print(f"    | J = {_j!r}")

    # Generate U-matrix
    u_mat_sh = []
    norb_sh = p['model']['norb_inequiv_sh']
    for norb, (_u, _up, _j) in zip(norb_sh, kanamori_sh):
        u_mat = numpy.zeros((norb, norb, norb, norb), numpy.complex_)
        for iorb, jorb in product(range(norb), repeat=2):
            u_mat[iorb, jorb, iorb, jorb] = _up
            u_mat[iorb, jorb, jorb, iorb] = _j
            u_mat[iorb, iorb, jorb, jorb] = _j
        for iorb in range(norb):
            u_mat[iorb, iorb, iorb, iorb] = _u
        u_mat_sh.append(u_mat)

    # Convert to spin-full U-matrix
    u_mat_so_sh = [to_spin_full_U_matrix(u_mat) for u_mat in u_mat_sh]

    return u_mat_so_sh


def _coefficients_ls_j(l, verbose=False, prefix=''):
    from sympy.physics.quantum.cg import CG
    from sympy import S

    s = S(1)/2

    # Lz, Sz
    m_s = numpy.array([[m, sz] for sz in [S(1)/2, -S(1)/2] for m in numpy.arange(-l, l+1)])

    # J, Jz
    # j_jz = numpy.array([[j, -j + _jz] for j in [l-s, l+s] for _jz in range(2*j+1)])
    j_jz = [[j, -j + _jz] for j in [l-s, l+s] for _jz in range((2*j+1)//2)]  # Jz<0
    j_jz = numpy.array(j_jz + [[j, -jz] for j, jz in j_jz])

    # matrix elements < L Lz S Sz | J Jz >
    mat_ls_j = numpy.zeros((len(j_jz), len(j_jz)), dtype=numpy.float)
    mat_for_print = []
    for i1, (m, ss) in enumerate(m_s):
        for i2, (j, jz) in enumerate(j_jz):
            cg = CG(l, m, s, ss, j, jz)
            me = cg.doit()
            # print(me)
            if me != 0:
                mat_ls_j[i1, i2] = me.evalf()
                mat_for_print.append((i1, i2, me))

    if verbose:
        print(f"{prefix}L = {l}")
        print(f"{prefix}S = {s}")
        print(f"{prefix}")
        print(f"{prefix}    Lz, Sz")
        for i, (m, sz) in enumerate(m_s):
            print(f"{prefix}{i:2d}: {m}, {sz}")
        print(f"{prefix}")
        print(f"{prefix}    J, Jz")
        for i, (j, jz) in enumerate(j_jz):
            print(f"{prefix}{i:2d}: {j}, {jz}")
        print(f"{prefix}")
        print(f"{prefix}(Lz, Sz)  (J, Jz)  < L Lz S Sz | J Jz >")
        for i1, i2, me in mat_for_print:
            print(f"{prefix}{i1:2d} {i2:2d}  {me}")

    # ls_basis = [f"l{m:+d},s{sz}" for m, sz in m_s]
    # j_basis = [f"j{j}{'+' if jz>0 else ''}{jz}" for j, jz in j_jz]
    # return mat_ls_j, ls_basis, j_basis

    return mat_ls_j, m_s, j_jz


def _from_ls_to_j(umat_ls, l, order=None):
    # print("\n Transform basis from LS to J")

    dim = 2*(2*l+1)
    assert umat_ls.shape == (dim, dim, dim, dim)

    # Get transformation matrix T
    tmat, basis_ls, basis_j = _coefficients_ls_j(l, verbose=False)
    # tmat, basis_ls, basis_j = _coefficients_ls_j(l, verbose=True, prefix='    ')

    assert tmat.shape == (dim, dim)
    assert basis_j.shape == (dim, 2)
    assert basis_ls.shape == (dim, 2)

    # Transform basis of U-matrix
    umat_j = numpy.einsum('mi,nj,ijkl,ko,lp', tmat.T.conj(), tmat.T.conj(), umat_ls, tmat, tmat)

    assert umat_j.shape == (dim, dim, dim, dim)

    return umat_j, basis_ls, basis_j


def _basis_names(l, basis):
    if basis == 'spherical':
        spdf = ['s', 'p', 'd', 'f']
        return numpy.array([f"{spdf[l]}{m:+d}" for m in range(-l, l+1)])
    elif basis == 'cubic':
        if l == 0:
            return numpy.array(('s',))
        else:
            return numpy.array(cubic_names(l))
    else:
        raise Exception("Here should not be arrived")


def _generate_umat_slater(p, l_sh, f_sh):
    #
    # Parse slater_basis
    #
    slater_basis = p['model']['slater_basis']
    nsh = p['model']['n_inequiv_shells']
    if slater_basis in ('cubic', 'spherical', 'spherical_j'):
        basis_sh = [slater_basis] * nsh
        order_sh = [None] * nsh
    else:
        slater_basis_sh = _parse_interaction_parameters(slater_basis, "slater_basis", nsh)
        basis_sh = [basis for basis, *_ in slater_basis_sh]
        order_sh = [order for _, *order in slater_basis_sh]
    print(f" slater_basis(basis) = {basis_sh!r}")
    print(f" slater_basis(order) = {order_sh!r}")

    # Check basis
    for basis in basis_sh:
        if basis not in ('cubic', 'spherical', 'spherical_j'):
            print(f"ERROR: basis={basis!r} not supported", file=sys.stderr)
            exit(1)

    # special treatment for basis='spherical_j'
    jbasis_sh = [False for _ in range(nsh)]
    for ish, basis in enumerate(basis_sh):
        if basis == 'spherical_j':
            basis_sh[ish] = 'spherical'  # replace
            jbasis_sh[ish] = True
    if numpy.any(numpy.array(jbasis_sh)) and p['model']['spin_orbit'] == False:
        print("Warning: 'spherical_j' in slater_basis should be used with spin_orbit=True", file=sys.stderr)

    # Support special symbols like order='eg', 't2g'
    for ish, (l, basis, order) in enumerate(zip(l_sh, basis_sh, order_sh)):
        if not order:  # None or []
            continue
        order_str = order[0]
        if isinstance(order_str, str):
            order = {
                (2, 'cubic', 'eg') : [2, 4],
                (2, 'cubic', 't2g') : [0, 1, 3],
            }.get((l, basis, order_str))
            if order is None:
                print(f"Error ! Unsupported pair of (l, basis, order) = ({l}, {basis}, {order_str})", file=sys.stderr)
                sys.exit(-1)
            order_sh[ish] = order  # replace
    # print(f"order_sh = {order_sh!r}")
    #
    # Generate U-matrix
    #
    u_mat_sh = []
    names_sh = []
    for l, f, basis, order in zip(l_sh, f_sh, basis_sh, order_sh):
        # basis names
        names = _basis_names(l=l, basis=basis)

        # U-matrix
        if l == 0:
            u_mat = numpy.full((1, 1, 1, 1), f[0], numpy.complex_)
        else:
            u_mat = U_matrix(l=l, radial_integrals=f, basis=basis)
        assert u_mat.shape == (len(names),) * 4

        u_mat_sh.append(u_mat)
        names_sh.append(names)

    # Change the order of bases
    # for ish, (u_mat, names, order, jbasis) in enumerate(zip(u_mat_sh, names_sh, order_sh, jbasis_sh)):
    #     if jbasis==False and order:  # exclude None, []
    #         u_mat_sh[ish] = u_mat[numpy.ix_(order, order, order, order)]
    #         names_sh[ish] = names[order]

    # print summary
    print("\n Slater interactions")
    for ish, (l, f, names) in enumerate(zip(l_sh, f_sh, names_sh)):
        print(f"  ish = {ish}")
        print(f"    | l = {l}")
        print(f"    | F_2m = {f}")
        print(f"    | basis = {names}")

    # Check the number of bases
    # norb_sh = p['model']['norb_inequiv_sh']
    # for ish, (names, norb) in enumerate(zip(names_sh, norb_sh)):
    #     if len(names) != norb:
    #         print(f"Error ! len(basis)={len(names)} is inconsistent with norb={norb} for ish={ish}")
    #         exit(1)
    #
    # Convert to spin-full U-matrix
    #
    u_mat_so_sh = [to_spin_full_U_matrix(u_mat) for u_mat in u_mat_sh]
    # print(names_sh)
    names_so_sh = [numpy.append(names, names) for names in names_sh]
    # print(names_so_sh)

    # Transform the basis from LS to J
    # names_so_sh = []
    for ish, (jbasis, u_mat_so, l) in enumerate(zip(jbasis_sh, u_mat_so_sh, l_sh)):
        if jbasis:
            # print(f"\n ish={ish} : Transform basis from LS to J")
            u_mat_so_sh[ish], _, basis_j = _from_ls_to_j(u_mat_so, l)
            # names_so_sh.append(names_so)
            # names_so_sh[ish] = names_so
            names_so_sh[ish] = numpy.array([f"j{j}{'+' if jz>0 else ''}{jz}" for j, jz in basis_j])  # convert to str

    # Change the order of bases
    for ish, (u_mat_so, names_so, order) in enumerate(zip(u_mat_so_sh, names_so_sh, order_sh)):
        # if jbasis and order:  # exclude None, []
        if order:  # exclude None, []
            dim = len(names_so)//2
            order_so = order + [i + dim for i in order]
            u_mat_so_sh[ish] = u_mat_so[numpy.ix_(order_so, order_so, order_so, order_so)]
            # print(names_so.shape)
            names_so_sh[ish] = names_so[order_so]

    # print summary
    # if numpy.any(numpy.array(jbasis_sh)):
    print("\n Basis in SO reps (after transformed, reordered, or trancated)")
    for ish, (jbasis, u_mat_so, names_so) in enumerate(zip(jbasis_sh, u_mat_so_sh, names_so_sh)):
        # names_j = [f"j{j}{'+' if jz>0 else ''}{jz}" for j, jz in names_so]
        print(f"  ish = {ish}")
        print(f"    | basis(up) = {names_so[:len(names_so)//2]}")
        print(f"    | basis(dn) = {names_so[len(names_so)//2:]}")

    # Check the number of bases
    norb_sh = p['model']['norb_inequiv_sh']
    for ish, (names_so, norb) in enumerate(zip(names_so_sh, norb_sh)):
        if len(names_so) != 2*norb:
            print(f"Error ! len(basis)={len(names_so)//2} is inconsistent with norb={norb} for ish={ish}")
            exit(1)

    return u_mat_so_sh


def _generate_umat_slater_uj(p):
    _check_parameters(p['model'], required=['slater_uj'], unused=['slater_f', 'kanamori'])

    def _U_J_to_F(_l, _u, _j):
        if _l == 0:
            return numpy.array([_u,], dtype=numpy.float_)
        else:
            return U_J_to_radial_integrals(_l, _u, _j)

    # parse kanamori parameters
    nsh = p['model']['n_inequiv_shells']
    slater_sh = _parse_interaction_parameters(p["model"]["slater_uj"], name='slater_uj', nsh=nsh, n_inner=3)
    l_sh = [int(l) for l, _, _ in slater_sh]
    f_sh = [_U_J_to_F(int(l), u, j) for l, u, j in slater_sh]

    # Generate U-matrix
    return _generate_umat_slater(p, l_sh, f_sh)


def _generate_umat_slater_f(p):
    _check_parameters(p['model'], required=['slater_f'], unused=['slater_uj', 'kanamori'])

    # parse slater parameters
    nsh = p['model']['n_inequiv_shells']
    slater_sh = _parse_interaction_parameters(p["model"]["slater_f"], name='slater_f', nsh=nsh, n_inner=5)
    l_sh = [int(l) for l, *_ in slater_sh]
    f_sh = [numpy.array(f[0:l+1], dtype=numpy.float_) for l, *f in slater_sh]

    # Warn if non-zero values are neglected from slater_f
    slater_f_neglected = [numpy.any(numpy.array(f[l+1:]) != 0) for l, *f in slater_sh]
    if numpy.any(numpy.array(slater_f_neglected)):
        print(f"Warning: Some non-zero values are neglected from slater_f={p['model']['slater_f']}. Only F_0, ..., F_2l are used.", file=sys.stderr)

    # Generate U-matrix
    return _generate_umat_slater(p, l_sh, f_sh)


def _generate_umat_respack(p):
    _check_parameters(p['model'], required=[], unused=['kanamori', 'slater_f', 'slater_uj'])

    nsh = p['model']['n_inequiv_shells']
    norb = p['model']['norb_inequiv_sh']
    #
    # Read U-matrix
    #
    w90u = Wannier90Converter(seedname=p["model"]["seedname"])
    w90u.convert_dft_input()
    nr_u, rvec_u, rdeg_u, nwan_u, hamr_u = w90u.read_wannier90hr(p["model"]["seedname"] + "_ur.dat")
    w90j = Wannier90Converter(seedname=p["model"]["seedname"])
    w90j.convert_dft_input()
    nr_j, rvec_j, rdeg_j, nwan_j, hamr_j = w90j.read_wannier90hr(p["model"]["seedname"] + "_jr.dat")
    #
    # Read 2-index U-matrix
    #
    umat2 = numpy.zeros((nwan_u, nwan_u), numpy.complex_)
    for ir in range(nr_u):
        if rvec_u[ir, 0] == 0 and rvec_u[ir, 1] == 0 and rvec_u[ir, 2] == 0:
            umat2 = hamr_u[ir]
    #
    # Read 2-index J-matrix
    #
    jmat2 = numpy.zeros((nwan_j, nwan_j), numpy.complex_)
    for ir in range(nr_j):
        if rvec_j[ir, 0] == 0 and rvec_j[ir, 1] == 0 and rvec_j[ir, 2] == 0:
            jmat2 = hamr_j[ir]
    #
    # Map into 4-index U at each correlated shell
    #
    u_mat_sh = [numpy.zeros((norb[ish], norb[ish], norb[ish], norb[ish]), numpy.complex_) for ish in range(nsh)]
    start = 0
    for ish in range(nsh):
        for iorb in range(norb[ish]):
            for jorb in range(norb[ish]):
                u_mat_sh[ish][iorb, jorb, iorb, jorb] = umat2[start+iorb, start+jorb]
                u_mat_sh[ish][iorb, jorb, jorb, iorb] = jmat2[start+iorb, start+jorb]
                u_mat_sh[ish][iorb, iorb, jorb, jorb] = jmat2[start+iorb, start+jorb]
        for iorb in range(norb[ish]):
            u_mat_sh[ish][iorb, iorb, iorb, iorb] = umat2[start+iorb, start+iorb]
        start += norb[ish]

    # Make U real
    # TODO: warn if U-matrix is not real
    for ish in range(nsh):
        u_mat_sh[ish][:, :, :, :].imag = 0.0

    # Convert to spin-full U-matrix
    u_mat_so_sh = [to_spin_full_U_matrix(u_mat) for u_mat in u_mat_sh]

    return u_mat_so_sh


def __generate_umat(p):
    """
    Add U-matrix block (Tentative)
    :param p: dictionary
        Input parameters
    :return:
    """
    # Add U-matrix block (Tentative)
    # ####  The format of this block is not fixed  ####
    #
    # Generate U-matrix
    #
    interaction = p["model"]["interaction"]
    func_umat = {
        'kanamori': _generate_umat_kanamori,
        'slater_uj': _generate_umat_slater_uj,
        'slater_f': _generate_umat_slater_f,
        'respack': _generate_umat_respack,
    }.get(interaction)
    if func_umat is None:
        print(f"Error ! Invalid interaction : {interaction}")
        sys.exit(-1)
    u_mat_so_sh = func_umat(p)
    #
    # Check U-matrix
    #
    nsh = p['model']['n_inequiv_shells']
    norb_sh = p['model']['norb_inequiv_sh']
    assert len(u_mat_so_sh) == nsh
    for u_mat, norb in zip(u_mat_so_sh, norb_sh):
        # assert u_mat.dtype == numpy.complex  # U-matrix is complex
        assert u_mat.shape == (2*norb, 2*norb, 2*norb, 2*norb)
    #
    # Convert to spin-full U-matrix
    #
    # u_mat_so_sh = [to_spin_full_U_matrix(u_mat) for u_mat in u_mat_sh]
    #
    # Transform LS basis to J
    #
    # if basis == 'spherical_j':
    if False:
        u_mat_so_sh = [_from_ls_to_j(u_mat_so, l) for u_mat_so, l in zip(u_mat_so_sh, l_sh)]
    #
    # Extract only density-density interactions if specified
    #
    if p["model"]["density_density"]:
        u_mat_so_sh = [umat2dd(u_mat_so) for u_mat_so in u_mat_so_sh]
    #
    # Write U-matrix
    #
    print("\n  @ Write the information of interactions")
    with HDFArchive(p["model"]["seedname"]+'.h5', 'a') as f:
        if "DCore" not in f:
            f.create_group("DCore")

        f["DCore"]["Umat"] = u_mat_so_sh
        print("\n    Written to {0}".format(p["model"]["seedname"]+'.h5'))


def __generate_local_potential(p):
    print("\n  @ Write the information of local potential")

    # str
    local_potential_matrix = p["model"]["local_potential_matrix"]
    local_potential_factor = p["model"]["local_potential_factor"]

    n_inequiv_shells = p["model"]['n_inequiv_shells']
    spin_orbit = p["model"]["spin_orbit"]

    # read parameters from DFT data
    skc = SumkDFTCompat(p["model"]["seedname"] + '.h5')

    assert skc.n_inequiv_shells == n_inequiv_shells

    corr_shells = skc.corr_shells
    dim_sh = [corr_shells[skc.inequiv_to_corr[ish]]['dim'] for ish in range(n_inequiv_shells)]

    # set factor
    try:
        fac = ast.literal_eval(local_potential_factor)
        if isinstance(fac, float) or isinstance(fac, int):
            fac = [float(fac)] * n_inequiv_shells
        elif isinstance(fac, list) or isinstance(fac, tuple):
            assert len(fac) == n_inequiv_shells
        else:
            raise Exception("local_potential_factor should be float or list of length %d" % n_inequiv_shells)
    except Exception as e:
        print("Error: local_potential_factor =", local_potential_factor)
        print(e)
        exit(1)

    # print factor
    print("fac =", fac)

    # set potential matrix
    pot = set_potential(local_potential_matrix, "local_potential_matrix", n_inequiv_shells, dim_sh, spin_orbit)

    for ish in range(n_inequiv_shells):
        pot[ish] *= fac[ish]

    # check if potential matrix is hermitian
    def is_hermitian(mat):
        return numpy.allclose(mat, mat.transpose().conj())
    try:
        for ish, pot_ish in enumerate(pot):
            for sp in range(pot_ish.shape[0]):
                assert is_hermitian(pot_ish[sp]), "potential matrix for ish={} sp={} is not hermitian".format(ish, sp)
    except AssertionError as e:
        print("Error:", e)
        exit(1)

    # write potential matrix
    with HDFArchive(p["model"]["seedname"] + '.h5', 'a') as f:
        f["DCore"]["LocalPotential"] = pot
    print("\n    Written to {0}".format(p["model"]["seedname"]+'.h5'))


def __check_if_Hk_is_hermite(h5file):
    with h5py.File(h5file, 'r') as f:
        hk = float_to_complex_array(f['/dft_input/hopping'][()])
        for ik in range(hk.shape[0]):
            for ib in range(hk.shape[1]):
                max_val = numpy.amax(numpy.abs(hk[ik,ib,:,:]))
                if max_val < 1e-8:
                    continue
                diff = numpy.amax(numpy.abs(hk[ik,ib,:,:] - hk[ik,ib,:,:].conjugate().transpose()))/max_val
                message = 'H(k) is not hermite at ik={} and iblock={}, relative diff is {}.' .format(ik, ib, diff)
                if diff > 1e-2:
                    raise RuntimeError('Error: {}'.format(message))
                elif diff > 1e-8:
                    print('Warning: {}'.format(message))


def dcore_pre(filename):
    """
    Main routine for the pre-processing tool

    Parameters
    ----------
    filename : string
        Input-file name
    """

    print("\n@@@@@@@@@@@@@@@@@@@  Reading Input File  @@@@@@@@@@@@@@@@@@@@\n")
    print("Input File Name : ", filename)
    #
    # Construct a parser with default values
    #
    pars = create_parser(['model'])
    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    #
    # Summary of input parameters
    #
    print("\n  @ Parameter summary")
    print("\n    [model] block")
    for k, v in list(p["model"].items()):
        print(f"      {k} = {v!r}")

    #
    # remove HDF5 file if exists
    #
    h5_file = p['model']['seedname'] + '.h5'
    if p['model']['lattice'] != 'external':
        if os.path.exists(h5_file):
            print("Removing the existing model HDF5 file...")
            os.remove(h5_file)

    #
    # Lattice information
    #   -> create h5_file/dft_input
    #
    print("\n@@@@@@@@@@@@@@@@@@@  Generate Model-HDF5 File  @@@@@@@@@@@@@@@@@@@@\n")
    lattice_model = create_lattice_model(p)
    lattice_model.generate_model_file()

    #
    # Interaction
    #   -> create h5_file/DCore/umat
    #
    print("\nGenerating U-matrix")
    __generate_umat(p)

    #
    # Local potential
    #   -> create h5_file/DCore/local_potential
    #
    print("\nGenerating local potential")
    __generate_local_potential(p)

    #
    # Check HDF5 file
    #
    print('')
    print('@@@@@@@@@@@@@@@@@@@ Check Model-HDF5 file @@@@@@@@@@@@@@@@@@@@')
    __check_if_Hk_is_hermite(h5_file)
    print_local_fields(h5_file)

    #
    # Finish
    #
    print("\n@@@@@@@@@@@@@@@@@@@@@@  Done  @@@@@@@@@@@@@@@@@@@@@@@@\n")

    raise_if_mpi_imported()


def run():
    from dcore.option_tables import generate_all_description
    import argparse
    from dcore.version import version, print_header

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_pre.py',
        description='pre script for dcore.',
        usage='$ dcore_pre input',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description()
    )
    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name."
                        )
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)
    dcore_pre(args.path_input_file)
