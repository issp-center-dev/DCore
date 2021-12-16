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
import re
import ast
import h5py
from h5 import HDFArchive
from dcore.program_options import create_parser
from triqs.operators.util.U_matrix import U_J_to_radial_integrals, U_matrix, eg_submatrix, t2g_submatrix

from dcore.converters.wannier90 import Wannier90Converter

from dcore.tools import *
from dcore.sumkdft_compat import SumkDFTCompat

from dcore.lattice_models import create_lattice_model
from dcore.lattice_models.tools import print_local_fields
from dcore.program_options import parse_parameters

def __print_paramter(p, param_name):
    print(param_name + " = " + str(p[param_name]))

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
    nsh = p['model']['n_inequiv_shells']
    kanamori = numpy.zeros((nsh, 3), numpy.float_)
    slater_f = numpy.zeros((nsh, 4), numpy.float_)
    slater_l = numpy.zeros(nsh, numpy.int_)
    #
    # Read interaction from input file
    #
    if p["model"]["interaction"] == 'kanamori':
        if p["model"]["kanamori"] == "None":
            print("Error ! Parameter \"kanamori\" is not specified.")
            sys.exit(-1)
        if p["model"]["slater_f"] != "None":
            print("Error ! Parameter \"slater_f\" is specified but is not used.")
            sys.exit(-1)
        if p["model"]["slater_uj"] != "None":
            print("Error ! Parameter \"slater_uj\" is specified but is not used.")
            sys.exit(-1)

        kanamori_list = re.findall(r'\(\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)', p["model"]["kanamori"])
        if len(kanamori_list) != nsh:
            print("\nError! The length of \"kanamori\" is wrong.")
            sys.exit(-1)

        try:
            for i, _list in enumerate(kanamori_list):
                _kanamori = [w for w in re.split(r'[)(,]', _list) if len(w) > 0]
                for j in range(3):
                    kanamori[i, j] = float(_kanamori[j])
        except RuntimeError:
            raise RuntimeError("Error ! Format of u_j is wrong.")
    elif p["model"]["interaction"] == 'slater_uj':
        if p["model"]["slater_uj"] == "None":
            print("Error ! Parameter \"slater_uj\" is not specified.")
            sys.exit(-1)
        if p["model"]["slater_f"] != "None":
            print("Error ! Parameter \"slater_f\" is specified but is not used.")
            sys.exit(-1)
        if p["model"]["kanamori"] != "None":
            print("Error ! Parameter \"kanamori\" is specified but is not used.")
            sys.exit(-1)

        f_list = re.findall(r'\(\s*\d+\s*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)',
                            p["model"]["slater_uj"])
        if len(f_list) != nsh:
            print("\nError! The length of \"slater_uj\" is wrong.")
            sys.exit(-1)

        try:
            for i, _list in enumerate(f_list):
                _slater = [w for w in re.split(r'[)(,]', _list) if len(w) > 0]
                slater_l[i] = int(_slater[0])
                slater_u = float(_slater[1])
                slater_j = float(_slater[2])
                if slater_l[i] == 0:
                    slater_f[i, 0] = slater_u
                else:
                    slater_f[i, 0:slater_l[i]+1] = U_J_to_radial_integrals(slater_l[i], slater_u, slater_j)
        except RuntimeError:
            raise RuntimeError("Error ! Format of u_j is wrong.")
    elif p["model"]["interaction"] == 'slater_f':
        if p["model"]["slater_f"] == "None":
            print("Error ! Parameter \"slater_f\" is not specified.")
            sys.exit(-1)
        if p["model"]["kanamori"] != "None":
            print("Error ! Parameter \"kanamori\" is specified but is not used.")
            sys.exit(-1)
        if p["model"]["slater_uj"] != "None":
            print("Error ! Parameter \"slater_uj\" is specified but is not used.")
            sys.exit(-1)

        f_list = re.findall(r'\(\s*\d+\s*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)',
                            p["model"]["slater_f"])
        slater_f = numpy.zeros((nsh, 4), numpy.float_)
        slater_l = numpy.zeros(nsh, numpy.int_)
        if len(f_list) != nsh:
            print("\nError! The length of \"slater_f\" is wrong.")
            sys.exit(-1)

        try:
            for i, _list in enumerate(f_list):
                _slater = [w for w in re.split(r'[)(,]', _list) if len(w) > 0]
                slater_l[i] = int(_slater[0])
                for j in range(4):
                    slater_f[i, j] = float(_slater[j])
        except RuntimeError:
            raise RuntimeError("Error ! Format of u_j is wrong.")
    elif p["model"]["interaction"] == 'respack':
        if p["model"]["kanamori"] != "None":
            print("Error ! Parameter \"kanamori\" is specified but is not used.")
            sys.exit(-1)
        if p["model"]["slater_uj"] != "None":
            print("Error ! Parameter \"slater_uj\" is specified but is not used.")
            sys.exit(-1)
        if p["model"]["slater_f"] != "None":
            print("Error ! Parameter \"slater_f\" is specified but is not used.")
            sys.exit(-1)
    else:
        print("Error ! Invalid interaction : ", p["model"]["interaction"])
        sys.exit(-1)

    #
    # Generate and Write U-matrix
    #
    print("\n  @ Write the information of interactions")
    f = HDFArchive(p["model"]["seedname"]+'.h5', 'a')

    norb = p['model']['norb_inequiv_sh']

    if not ("DCore" in f):
        f.create_group("DCore")
    #
    u_mat = [numpy.zeros((norb[ish], norb[ish], norb[ish], norb[ish]), numpy.complex_) for ish in range(nsh)]
    if p["model"]["interaction"] == 'kanamori':
        for ish in range(nsh):
            for iorb in range(norb[ish]):
                for jorb in range(norb[ish]):
                    u_mat[ish][iorb, jorb, iorb, jorb] = kanamori[ish, 1]
                    u_mat[ish][iorb, jorb, jorb, iorb] = kanamori[ish, 2]
                    u_mat[ish][iorb, iorb, jorb, jorb] = kanamori[ish, 2]
            for iorb in range(norb[ish]):
                u_mat[ish][iorb, iorb, iorb, iorb] = kanamori[ish, 0]
    elif p["model"]["interaction"] == 'slater_uj' or p["model"]["interaction"] == 'slater_f':
        for ish in range(nsh):
            if slater_l[ish] == 0:
                umat_full = numpy.zeros((1, 1, 1, 1), numpy.complex_)
                umat_full[0, 0, 0, 0] = slater_f[ish, 0]
            else:
                umat_full = U_matrix(l=slater_l[ish], radial_integrals=slater_f[ish, :], basis='cubic')
            #
            # For t2g or eg, compute submatrix
            #
            if slater_l[ish]*2+1 != norb[ish]:
                if slater_l[ish] == 2 and norb[ish] == 2:
                    u_mat[ish] = eg_submatrix(umat_full)
                elif slater_l[ish] == 2 and norb[ish] == 3:
                    u_mat[ish] = t2g_submatrix(umat_full)
                else:
                    print("Error ! Unsupported pair of l and norb : ", slater_l[ish], norb[ish])
                    sys.exit(-1)
            else:
                u_mat[ish] = umat_full
    elif p["model"]["interaction"] == 'respack':
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
    for ish in range(nsh):
        u_mat[ish][:, :, :, :].imag = 0.0
    #
    # Spin & Orb
    #
    u_mat2 = [numpy.zeros((norb[ish]*2, norb[ish]*2, norb[ish]*2, norb[ish]*2), numpy.complex_)
              for ish in range(nsh)]
    for ish in range(nsh):
        u_mat2[ish] = to_spin_full_U_matrix(u_mat[ish])

    if p["model"]["density_density"]:
        for ish in range(nsh):
            umat = umat2dd(u_mat2[ish][:])
            u_mat2[ish][:] = umat
    f["DCore"]["Umat"] = u_mat2
    print("\n    Written to {0}".format(p["model"]["seedname"]+'.h5'))
    del f


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
        print("      {0} = {1}".format(k, v))

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
    __generate_umat(p)

    #
    # Local potential
    #   -> create h5_file/DCore/local_potential
    #
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
