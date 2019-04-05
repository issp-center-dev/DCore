#!/usr/bin/env python
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
from __future__ import print_function
import sys
import os
import numpy
import re
import ast
from pytriqs.archive.hdf_archive import HDFArchive
from program_options import create_parser
from pytriqs.operators.util.U_matrix import U_J_to_radial_integrals, U_matrix, eg_submatrix, t2g_submatrix

from converters.wannier90_converter import Wannier90Converter

from .tools import *

from .lattice_model import create_lattice_model, print_local_fields

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
    ncor = p["model"]['ncor']
    kanamori = numpy.zeros((ncor, 3), numpy.float_)
    slater_f = numpy.zeros((ncor, 4), numpy.float_)
    slater_l = numpy.zeros(ncor, numpy.int_)
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
        if len(kanamori_list) != ncor:
            print("\nError! The length of \"kanamori\" is wrong.")
            sys.exit(-1)

        try:
            for i, _list in enumerate(kanamori_list):
                _kanamori = filter(lambda w: len(w) > 0, re.split(r'[)(,]', _list))
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
        if len(f_list) != ncor:
            print("\nError! The length of \"slater_uj\" is wrong.")
            sys.exit(-1)

        try:
            for i, _list in enumerate(f_list):
                _slater = filter(lambda w: len(w) > 0, re.split(r'[)(,]', _list))
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
        slater_f = numpy.zeros((ncor, 4), numpy.float_)
        slater_l = numpy.zeros(ncor, numpy.int_)
        if len(f_list) != ncor:
            print("\nError! The length of \"slater_f\" is wrong.")
            sys.exit(-1)

        try:
            for i, _list in enumerate(f_list):
                _slater = filter(lambda w: len(w) > 0, re.split(r'[)(,]', _list))
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

    corr_shells = f["dft_input"]["corr_shells"]
    norb = [corr_shells[icor]["dim"] for icor in range(ncor)]
    if p["model"]["spin_orbit"] or p["model"]["non_colinear"]:
        for icor in range(ncor):
            norb[icor] //= 2

    if not ("DCore" in f):
        f.create_group("DCore")
    #
    u_mat = [numpy.zeros((norb[icor], norb[icor], norb[icor], norb[icor]), numpy.complex_) for icor in range(ncor)]
    if p["model"]["interaction"] == 'kanamori':
        for icor in range(ncor):
            for iorb in range(norb[icor]):
                for jorb in range(norb[icor]):
                    u_mat[icor][iorb, jorb, iorb, jorb] = kanamori[icor, 1]
                    u_mat[icor][iorb, jorb, jorb, iorb] = kanamori[icor, 2]
                    u_mat[icor][iorb, iorb, jorb, jorb] = kanamori[icor, 2]
            for iorb in range(norb[icor]):
                u_mat[icor][iorb, iorb, iorb, iorb] = kanamori[icor, 0]
    elif p["model"]["interaction"] == 'slater_uj' or p["model"]["interaction"] == 'slater_f':
        for icor in range(ncor):
            if slater_l[icor] == 0:
                umat_full = numpy.zeros((1, 1, 1, 1), numpy.complex_)
                umat_full[0, 0, 0, 0] = slater_f[icor, 0]
            else:
                umat_full = U_matrix(l=slater_l[icor], radial_integrals=slater_f[icor, :], basis='cubic')
            #
            # For t2g or eg, compute submatrix
            #
            if slater_l[icor]*2+1 != norb[icor]:
                if slater_l[icor] == 2 and norb[icor] == 2:
                    u_mat[icor] = eg_submatrix(umat_full)
                elif slater_l[icor] == 2 and norb[icor] == 3:
                    u_mat[icor] = t2g_submatrix(umat_full)
                else:
                    print("Error ! Unsupported pair of l and norb : ", slater_l[icor], norb[icor])
                    sys.exit(-1)
            else:
                u_mat[icor] = umat_full
    elif p["model"]["interaction"] == 'respack':
        #w90u = converters.wannier90_converter.Wannier90Converter(seedname=p["model"]["seedname"])
        w90u = Wannier90Converter(seedname=p["model"]["seedname"])
        nr_u, rvec_u, rdeg_u, nwan_u, hamr_u = w90u.read_wannier90hr(p["model"]["seedname"] + "_ur.dat")
        #w90j = converters.wannier90_converter.Wannier90Converter(seedname=p["model"]["seedname"])
        w90j = Wannier90Converter(seedname=p["model"]["seedname"])
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
        for icor in range(ncor):
            for iorb in range(norb[icor]):
                for jorb in range(norb[icor]):
                    u_mat[icor][iorb, jorb, iorb, jorb] = umat2[start+iorb, start+jorb]
                    u_mat[icor][iorb, jorb, jorb, iorb] = jmat2[start+iorb, start+jorb]
                    u_mat[icor][iorb, iorb, jorb, jorb] = jmat2[start+iorb, start+jorb]
            for iorb in range(norb[icor]):
                u_mat[icor][iorb, iorb, iorb, iorb] = umat2[start+iorb, start+iorb]
            start += norb[icor]
    #
    for icor in range(ncor):
        u_mat[icor][:, :, :, :].imag = 0.0
    #
    # Spin & Orb
    #
    u_mat2 = [numpy.zeros((norb[icor]*2, norb[icor]*2, norb[icor]*2, norb[icor]*2), numpy.complex_)
              for icor in range(ncor)]
    for icor in range(ncor):
        u_mat2[icor] = to_spin_full_U_matrix(u_mat[icor])
    f["DCore"]["Umat"] = u_mat2
    print("\n    Wrote to {0}".format(p["model"]["seedname"]+'.h5'))
    del f


def __generate_local_potential(p):
    print("\n  @ Write the information of local potential")

    # str
    local_potential_matrix = p["model"]["local_potential_matrix"]
    local_potential_factor = p["model"]["local_potential_factor"]

    ncor = p["model"]['ncor']
    spin_orbit = p["model"]["spin_orbit"]

    # read parameters from DFT data
    with HDFArchive(p["model"]["seedname"] + '.h5', 'r') as f:
        corr_shells = f["dft_input"]["corr_shells"]
    dim_crsh = [corr_shell["dim"] for corr_shell in corr_shells]

    # set factor
    try:
        fac = ast.literal_eval(local_potential_factor)
        if isinstance(fac, float):
            fac = [fac] * ncor
        elif isinstance(fac, list) or isinstance(fac, tuple):
            assert len(fac) == ncor
        else:
            raise Exception("local_potential_factor should be float or list of length %d" % ncor)
    except Exception as e:
        print("Error: local_potential_factor =", local_potential_factor)
        print(e)
        exit(1)

    # print factor
    print("fac =", fac)

    # init potential
    # pot.shape = (2, orb1, orb2)     w/  spin-orbit
    #             (1, 2*orb1, 2*orb2) w/o spin-orbit
    if not spin_orbit:
        pot = [numpy.zeros((2, dim_crsh[icor], dim_crsh[icor]), numpy.complex_) for icor in range(ncor)]
    else:
        pot = [numpy.zeros((1, dim_crsh[icor], dim_crsh[icor]), numpy.complex_) for icor in range(ncor)]

    # read potential matrix
    if local_potential_matrix != 'None':
        try:
            files = ast.literal_eval(local_potential_matrix)
            assert isinstance(files, dict)
            assert all([ish < ncor for ish in files.keys()])
        except Exception as e:
            print("Error: local_potential_matrix =", local_potential_matrix)
            print(e)
            exit(1)

        for ish, file in files.items():
            read_potential(file, pot[ish])
            pot[ish] *= fac[ish]

    # print potential
    print("\n--- potential matrix")
    for ish, pot_ish in enumerate(pot):
        print("ish =", ish)
        for sp in range(pot_ish.shape[0]):
            print("sp =", sp)
            print(pot_ish[sp])

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
    print("\n    Wrote to {0}".format(p["model"]["seedname"]+'.h5'))


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
    pars = create_parser()
    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    ncor = p["model"]['ncor']
    #
    # Summary of input parameters
    #
    print("\n  @ Parameter summary")
    print("\n    [model] block")
    for k, v in p["model"].items():
        print("      {0} = {1}".format(k, v))
    print("\n    [system] block")
    for k, v in p["system"].items():
        print("      {0} = {1}".format(k, v))

    if os.path.exists(p['model']['seedname'] + '.h5'):
        print("Removing the existing model HDF5 file...")
        os.remove(p['model']['seedname'] + '.h5')

    #
    # One-body term
    #
    print("\n@@@@@@@@@@@@@@@@@@@  Generate Model-HDF5 File  @@@@@@@@@@@@@@@@@@@@\n")
    lattice_model = create_lattice_model(p)
    lattice_model.generate_model_file()

    # Interaction
    #
    __generate_umat(p)

    #
    # Local potential
    #
    __generate_local_potential(p)

    #
    # Check one-body term
    #
    print('')
    print('@@@@@@@@@@@@@@@@@@@ Check Model-HDF5 file @@@@@@@@@@@@@@@@@@@@')
    print_local_fields(p['model']['seedname'] + '.h5')

    #
    # Finish
    #

    print("\n@@@@@@@@@@@@@@@@@@@@@@  Done  @@@@@@@@@@@@@@@@@@@@@@@@\n")


if __name__ == '__main__':
    from .option_tables import generate_all_description
    import argparse

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

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file is not exist.")
        sys.exit(-1)
    dcore_pre(args.path_input_file)
