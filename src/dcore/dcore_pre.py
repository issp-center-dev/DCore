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

from dcore._dispatcher import HDFArchive
from dcore.program_options import create_parser

from dcore.tools import *
from dcore.sumkdft_compat import SumkDFTCompat

from dcore.lattice_models import create_lattice_model
from dcore.lattice_models.tools import print_local_fields
from dcore.program_options import parse_parameters, print_parameters, delete_parameters
from dcore.interaction import generate_umat


def __generate_local_potential(p):
    print("\n  @ Write the information of local potential")

    # str
    local_potential_matrix = p["model"]["local_potential_matrix"]
    local_potential_factor = p["model"]["local_potential_factor"]

    n_inequiv_shells = p["model"]['n_inequiv_shells']
    spin_orbit = p["model"]["spin_orbit"]

    # read parameters from DFT data
    try:
        skc = SumkDFTCompat(p["model"]["seedname"] + '.h5')
    except RuntimeError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        print("Generate lattice model by running dcore_pre with skip_lattice=False (default) before generating local potential.", file=sys.stderr)
        sys.exit(1)

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
        print("Error: local_potential_factor =", local_potential_factor, file=sys.stderr)
        print(e, file=sys.stderr)
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
        print("Error:", e, file=sys.stderr)
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


def dcore_pre(input_filenames):
    """
    Main routine for the pre-processing tool

    Parameters
    ----------
    filename : string
        Input-file name
    """

    print("\n@@@@@@@@@@@@@@@@@@@  Reading Input File  @@@@@@@@@@@@@@@@@@@@\n")
    print("Input Filename(s) : ", input_filenames)
    #
    # Construct a parser with default values
    #
    pars = create_parser(['model', 'pre'])
    #
    # Parse keywords and store
    #
    pars.read(input_filenames)
    p = pars.as_dict()
    parse_parameters(p)

    #
    # Delete unnecessary parameters
    #
    delete_parameters(p, block='model', delete=['bvec'])

    #
    # Summary of input parameters
    #
    print_parameters(p)

    print("\n@@@@@@@@@@@@@@@@@@@  Generate Model-HDF5 File  @@@@@@@@@@@@@@@@@@@@")
    #
    # remove HDF5 file if exists
    #
    h5_file = p['model']['seedname'] + '.h5'
    if os.path.exists(h5_file):
        print(f"\nFile '{h5_file}' found")
        if p['pre']['overwrite']:
            print("  --> overwritten (overwrite=True)")
        else:
            print("  --> replaced (overwrite=False)")
            os.remove(h5_file)
    else:
        print(f"\nFile '{h5_file}' is created")


    #
    # Lattice information
    #   -> create h5_file/dft_input
    #
    print("\nGenerating lattice model including H(k)")
    print(f"  in {h5_file}/dft_input")
    if p['pre']['skip_lattice']:
        print("skip")
    else:
        lattice_model = create_lattice_model(p)
        lattice_model.generate_model_file()

    #
    # Interaction
    #   -> create h5_file/DCore/umat
    #
    print("\nGenerating U-matrix")
    print(f"  in {h5_file}/DCore/Umat")
    if p['pre']['skip_umat']:
        print("skip")
    else:
        generate_umat(p)

    #
    # Local potential
    #   -> create h5_file/DCore/local_potential
    #
    print("\nGenerating local potential")
    print(f"  in {h5_file}/DCore/LocalPotential")
    if p['pre']['skip_local_potential']:
        print("skip")
    else:
        __generate_local_potential(p)

    #
    # Check HDF5 file
    #
    print('')
    print('@@@@@@@@@@@@@@@@@@@ Check Model-HDF5 file @@@@@@@@@@@@@@@@@@@@')

    print("\nChecking H(k)")
    if p['pre']['skip_lattice']:
        print("skip")
    else:
        __check_if_Hk_is_hermite(h5_file)

    print("\nLocal Fields")
    if p['pre']['skip_lattice']:
        print("skip")
    else:
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
        description='Pre-processing script in DCore',
        # usage='$ dcore_pre input',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description()
    )
    parser.add_argument('path_input_files',
                        action='store',
                        default=None,
                        type=str,
                        nargs='*',
                        help="Input filename(s)",
                        )
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))

    args = parser.parse_args()
    for path_input_file in args.path_input_files:
        if os.path.isfile(path_input_file) is False:
            sys.exit(f"Input file '{path_input_file}' does not exist.")
    dcore_pre(args.path_input_files)
