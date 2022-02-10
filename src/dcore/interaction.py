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
import numpy
import ast
from typing import List, Dict, Optional, Union, Sequence, cast

from itertools import product
from dcore._dispatcher import HDFArchive, U_J_to_radial_integrals, U_matrix, cubic_names

from dcore.converters.wannier90 import Wannier90Converter

from dcore.tools import *


def _check_parameters(p: Dict, required: List[str], unused: List[str]):
    assert isinstance(p, dict)
    assert isinstance(required, list)
    assert isinstance(unused, list)
    for key in required:
        if p[key] == "None":
            print(f"Error ! Parameter '{key}' is not specified.", file=sys.stderr)
            sys.exit(-1)
    for key in unused:
        if p[key] != "None":
            print(f"Error ! Parameter '{key}' is specified but is not used.", file=sys.stderr)
            sys.exit(-1)


def _parse_interaction_parameters(
    input_str: str, name, nsh: int, n_inner: Optional[int]=None) -> Sequence[Sequence[Union[int,float,str]]]:
    # parse interaction parameters
    # return list of list
    assert isinstance(input_str, str)
    assert isinstance(name, str)
    assert isinstance(nsh, int)
    assert n_inner is None or isinstance(n_inner, int)
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


def _generate_umat_kanamori(p: Dict):
    assert isinstance(p, Dict)
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
        u_mat = numpy.zeros((norb, norb, norb, norb), numpy.complex128)
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


def _coefficients_ls_j(l: int, verbose: bool=False, prefix: str=''):
    assert isinstance(l, int)
    assert isinstance(verbose, bool)
    assert isinstance(prefix, str)
    from sympy.physics.quantum.cg import CG
    from sympy import S

    # Lz, Sz
    s = S(1)/2
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
            coef = cg.doit()
            # print(coef)
            if coef != 0:
                mat_ls_j[i1, i2] = coef.evalf()
                mat_for_print.append((i1, i2, coef))

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
        for i1, i2, coef in mat_for_print:
            print(f"{prefix}{i1:2d} {i2:2d}  {coef}")

    return mat_ls_j, m_s, j_jz


def _from_ls_to_j(umat_ls: numpy.ndarray, l: int, order=None):
    assert isinstance(umat_ls, numpy.ndarray)
    assert isinstance(l, int)
    dim = 2*(2*l+1)
    assert umat_ls.shape == (dim, dim, dim, dim)

    # Get transformation matrix T
    tmat, basis_ls, basis_j = _coefficients_ls_j(l, verbose=False)

    assert tmat.shape == (dim, dim)
    assert basis_j.shape == (dim, 2)
    assert basis_ls.shape == (dim, 2)

    # Transform basis of U-matrix
    umat_j = numpy.einsum('mi,nj,ijkl,ko,lp', tmat.T.conj(), tmat.T.conj(), umat_ls, tmat, tmat)

    assert umat_j.shape == (dim, dim, dim, dim)

    return umat_j, basis_ls, basis_j


def _basis_names(l: int, basis: str):
    assert isinstance(l, int)
    assert isinstance(basis, str)
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


def _generate_umat_slater(p: Dict, l_sh: List[int], f_sh: List[numpy.ndarray]):
    #
    # Parse slater_basis
    #
    assert isinstance(l_sh, list) and all((isinstance(l_, int) for l_ in l_sh))
    assert isinstance(f_sh, list) and all((isinstance(f_, numpy.ndarray) for f_ in f_sh))

    slater_basis = p['model']['slater_basis']
    nsh = p['model']['n_inequiv_shells']
    if slater_basis in ('cubic', 'spherical', 'spherical_j'):
        basis_sh = [slater_basis] * nsh #type: List[str]
        order_sh = [None] * nsh #type: List[Optional[List]]
    else:
        slater_basis_sh = _parse_interaction_parameters(slater_basis, "slater_basis", nsh)
        basis_sh = [cast(str, basis) for basis, *_ in slater_basis_sh] #type: List[str]
        order_sh = [order for _, *order in slater_basis_sh] #type: List[Optional[List]]
    print(f" slater_basis(basis) = {basis_sh!r}")
    print(f" slater_basis(order) = {order_sh!r}")

    if 'spherical_j' in basis_sh and p['model']['spin_orbit'] == False:
        print("Warning: 'spherical_j' in slater_basis should be used with spin_orbit=True", file=sys.stderr)
    #
    # Generate U-matrix
    #
    u_mat_so_sh = []  # to be returned
    norb_sh = p['model']['norb_inequiv_sh']
    print("\n Slater interactions")
    for ish, (l, f, basis, order, norb) in enumerate(zip(l_sh, f_sh, basis_sh, order_sh, norb_sh)):
        print(f"  ish = {ish}")
        print(f"    | l = {l}")
        print(f"    | F_2m = {f}")

        # Check basis
        if basis not in ('cubic', 'spherical', 'spherical_j'):
            print(f"ERROR: basis={basis!r} not supported", file=sys.stderr)
            exit(1)

        # special treatment for basis='spherical_j'
        jbasis = False
        if basis == 'spherical_j':
            basis = 'spherical'  # replace
            jbasis = True

        # Support special symbols like order='eg', 't2g'
        if order:  # exclude None, []
            order_str = order[0]
            if isinstance(order_str, str):
                order = {
                    (2, 'cubic', 'eg') : [2, 4],
                    (2, 'cubic', 't2g') : [0, 1, 3],
                }.get((l, basis, order_str))
                if order is None:
                    print(f"Error ! Unsupported pair of (l, basis, order) = ({l!r}, {basis!r}, {order_str!r})", file=sys.stderr)
                    sys.exit(-1)
            _norb = len(order)
        else:
            _norb = 2*l + 1

        # Check the number of bases
        if norb != _norb:
            print(f"Error ! norb={norb} is inconsistent with (# of basis/sp)={_norb}", file=sys.stderr)
            exit(1)

        # Generate U-matrix
        if l == 0:
            u_mat = numpy.full((1, 1, 1, 1), f[0], numpy.complex128)
        else:
            u_mat = U_matrix(l=l, radial_integrals=f, basis=basis)

        # basis names
        names = _basis_names(l=l, basis=basis)
        assert u_mat.shape == (len(names),) * 4
        print(f"    | basis/sp = {names}")

        # Convert to spin-full U-matrix
        u_mat_so = to_spin_full_U_matrix(u_mat)
        names_so = numpy.append(names, names)

        # Transform the basis from LS to J
        if jbasis:
            u_mat_so, _, basis_j = _from_ls_to_j(u_mat_so, l)
            names_so = numpy.array([f"j{j}{'+' if jz>0 else ''}{jz}" for j, jz in basis_j])  # convert to str

        # Change the order of bases
        if order:  # exclude None, []
            order_so = order + [i + (2*l + 1) for i in order]
            u_mat_so = u_mat_so[numpy.ix_(order_so, order_so, order_so, order_so)]
            names_so = names_so[order_so]

        # Print summary in SO rep
        print(f"    |")
        print(f"    | in SO rep (after transformed, reordered, or truncated)")
        print(f"    | basis(up) = {names_so[:len(names_so)//2]}")
        print(f"    | basis(dn) = {names_so[len(names_so)//2:]}")

        assert len(names_so) == 2*norb
        assert u_mat_so.shape == (2*norb,) * 4

        u_mat_so_sh.append(u_mat_so)

    return u_mat_so_sh


def _generate_umat_slater_uj(p: Dict):
    assert isinstance(p, Dict)
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


def _generate_umat_slater_f(p: Dict):
    assert isinstance(p, Dict)
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


def _generate_umat_respack(p: Dict):
    assert isinstance(p, Dict)
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
    umat2 = numpy.zeros((nwan_u, nwan_u), numpy.complex128)
    for ir in range(nr_u):
        if rvec_u[ir, 0] == 0 and rvec_u[ir, 1] == 0 and rvec_u[ir, 2] == 0:
            umat2 = hamr_u[ir]
    #
    # Read 2-index J-matrix
    #
    jmat2 = numpy.zeros((nwan_j, nwan_j), numpy.complex128)
    for ir in range(nr_j):
        if rvec_j[ir, 0] == 0 and rvec_j[ir, 1] == 0 and rvec_j[ir, 2] == 0:
            jmat2 = hamr_j[ir]
    #
    # Map into 4-index U at each correlated shell
    #
    u_mat_sh = [numpy.zeros((norb[ish], norb[ish], norb[ish], norb[ish]), numpy.complex128) for ish in range(nsh)]
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


def generate_umat(p: Dict):
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
    assert isinstance(p, Dict)
    interaction = p["model"]["interaction"]
    func_umat = {
        'kanamori': _generate_umat_kanamori,
        'slater_uj': _generate_umat_slater_uj,
        'slater_f': _generate_umat_slater_f,
        'respack': _generate_umat_respack,
    }.get(interaction)
    if func_umat is None:
        print(f"Error ! Invalid interaction : {interaction}", file=sys.stderr)
        sys.exit(-1)
    u_mat_so_sh = func_umat(p)
    #
    # Check U-matrix
    #
    nsh = p['model']['n_inequiv_shells']
    norb_sh = p['model']['norb_inequiv_sh']
    assert len(u_mat_so_sh) == nsh
    for u_mat, norb in zip(u_mat_so_sh, norb_sh):
        # TODO: what to do if U-matrix is real
        # assert u_mat.dtype == numpy.complex  # U-matrix is complex
        assert u_mat.shape == (2*norb, 2*norb, 2*norb, 2*norb)
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