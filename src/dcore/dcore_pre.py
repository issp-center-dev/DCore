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
from dcore._dispatcher import HDFArchive, U_J_to_radial_integrals, U_matrix, eg_submatrix, t2g_submatrix
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


def _parse_interaction_parameters(input_str, name, nsh, n_inner):
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

    # print Kanamori parameters
    print("\n Kanamori parameters")
    print("  ish :  U,  U',  J")
    for ish, (u, up, j) in enumerate(kanamori_sh):
        print(f"    {ish} : {u!r}, {up!r}, {j!r}")

    # Generate U-matrix
    norb = p['model']['norb_inequiv_sh']
    u_mat = [numpy.zeros((norb[ish], norb[ish], norb[ish], norb[ish]), numpy.complex_) for ish in range(nsh)]
    for ish in range(nsh):
        for iorb, jorb in product(range(norb[ish]), repeat=2):
            u_mat[ish][iorb, jorb, iorb, jorb] = kanamori_sh[ish][1]
            u_mat[ish][iorb, jorb, jorb, iorb] = kanamori_sh[ish][2]
            u_mat[ish][iorb, iorb, jorb, jorb] = kanamori_sh[ish][2]
        for iorb in range(norb[ish]):
            u_mat[ish][iorb, iorb, iorb, iorb] = kanamori_sh[ish][0]
    return u_mat


def _generate_umat_slater(slater_l, slater_f, nsh, norb):
    # print F values
    print("\n Slater integrals F")
    print("  ish : l,  [F0 F2 F4 F6]")
    for ish in range(nsh):
        print(f"    {ish} : {slater_l[ish]},  {slater_f[ish]}")

    # Generate U-matrix
    u_mat = []
    for ish in range(nsh):
        if slater_l[ish] == 0:
            umat_full = numpy.zeros((1, 1, 1, 1), numpy.complex_)
            umat_full[0, 0, 0, 0] = slater_f[ish][0]
        else:
            umat_full = U_matrix(l=slater_l[ish], radial_integrals=slater_f[ish], basis='cubic')
        #
        # For t2g or eg, compute submatrix
        #
        if slater_l[ish]*2+1 != norb[ish]:
            if slater_l[ish] == 2 and norb[ish] == 2:
                umat_sub = eg_submatrix(umat_full)
            elif slater_l[ish] == 2 and norb[ish] == 3:
                umat_sub = t2g_submatrix(umat_full)
            else:
                print("Error ! Unsupported pair of l and norb : ", slater_l[ish], norb[ish])
                sys.exit(-1)
            u_mat.append(umat_sub)
        else:
            u_mat.append(umat_full)
    return u_mat


def _generate_umat_slater_uj(p):
    _check_parameters(p['model'], required=['slater_uj'], unused=['slater_f', 'kanamori'])

    def _U_J_to_F(_l, _u, _j):
        _slater_f = numpy.zeros(4, numpy.float_)
        if _l == 0:
            _slater_f[0] = _u
        else:
            _slater_f[0:_l+1] = U_J_to_radial_integrals(_l, _u, _j)
        return _slater_f

    # parse kanamori parameters
    nsh = p['model']['n_inequiv_shells']
    slater_sh = _parse_interaction_parameters(p["model"]["slater_uj"], name='slater_uj', nsh=nsh, n_inner=3)
    slater_l = [int(l) for l, _, _ in slater_sh]
    slater_f = [_U_J_to_F(int(l), u, j) for l, u, j in slater_sh]

    # Generate U-matrix
    norb = p['model']['norb_inequiv_sh']
    return _generate_umat_slater(slater_l, slater_f, nsh, norb)


def _generate_umat_slater_f(p):
    _check_parameters(p['model'], required=['slater_f'], unused=['slater_uj', 'kanamori'])

    # parse kanamori parameters
    nsh = p['model']['n_inequiv_shells']
    slater_sh = _parse_interaction_parameters(p["model"]["slater_f"], name='slater_f', nsh=nsh, n_inner=5)
    slater_l = [int(slater[0]) for slater in slater_sh]
    slater_f = [numpy.array(slater[1:], dtype=numpy.float_) for slater in slater_sh]

    # Generate U-matrix
    norb = p['model']['norb_inequiv_sh']
    return _generate_umat_slater(slater_l, slater_f, nsh, norb)


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
    u_mat = [numpy.zeros((norb[ish], norb[ish], norb[ish], norb[ish]), numpy.complex_) for ish in range(nsh)]
    start = 0
    for ish in range(nsh):
        for iorb in range(norb[ish]):
            for jorb in range(norb[ish]):
                u_mat[ish][iorb, jorb, iorb, jorb] = umat2[start+iorb, start+jorb]
                u_mat[ish][iorb, jorb, jorb, iorb] = jmat2[start+iorb, start+jorb]
                u_mat[ish][iorb, iorb, jorb, jorb] = jmat2[start+iorb, start+jorb]
        for iorb in range(norb[ish]):
            u_mat[ish][iorb, iorb, iorb, iorb] = umat2[start+iorb, start+iorb]
        start += norb[ish]

    # Make U real
    # TODO: warn if U-matrix is not real
    for ish in range(nsh):
        u_mat[ish][:, :, :, :].imag = 0.0

    return u_mat


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
    if interaction == 'kanamori':
        u_mat = _generate_umat_kanamori(p)
    elif interaction == 'slater_uj':
        u_mat = _generate_umat_slater_uj(p)
    elif interaction == 'slater_f':
        u_mat = _generate_umat_slater_f(p)
    elif interaction == 'respack':
        u_mat = _generate_umat_respack(p)
    else:
        print(f"Error ! Invalid interaction : {interaction}")
        sys.exit(-1)
    #
    # Check U-matrix
    #
    nsh = p['model']['n_inequiv_shells']
    assert len(u_mat) == nsh
    for ish in range(nsh):
        assert u_mat[ish].dtype == numpy.complex  # U-matrix is complex
    #
    # Convert to spin-full U-matrix
    #
    u_mat2 = [to_spin_full_U_matrix(u_mat[ish]) for ish in range(nsh)]
    #
    # Extract only density-density interactions if specified
    #
    if p["model"]["density_density"]:
        for ish in range(nsh):
            umat = umat2dd(u_mat2[ish][:])
            u_mat2[ish][:] = umat
    #
    # Write U-matrix
    #
    print("\n  @ Write the information of interactions")
    with HDFArchive(p["model"]["seedname"]+'.h5', 'a') as f:
        if "DCore" not in f:
            f.create_group("DCore")

        f["DCore"]["Umat"] = u_mat2
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
