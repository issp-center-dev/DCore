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


from .typed_parser import *

import numpy
import re
import ast
from collections import namedtuple

from .mpi_command import default_mpi_command


def create_parser(target_sections=None):
    """
    Create a parser for all program options of DCore
    """
    if target_sections is None:
        parser = TypedParser(['mpi', 'model', 'system', 'impurity_solver', 'control', 'tool', 'bse'])
    else:
        parser = TypedParser(target_sections)

    # [mpi]
    parser.add_option("mpi", "command", str, default_mpi_command, "Command for executing a MPI job. # will be relaced by the number of processes.")

    # [model]
    parser.add_option("model", "seedname", str, "dcore", "Name of the system. The model HDF5 file will be seedname.h5.")
    parser.add_option("model", "lattice", str, "chain",
                      'Chosen from "chain", "square", "cubic", "bethe", "wannier90", and "external"')
    parser.add_option("model", "t", float, 1.0, "Transfer integral (Nearest neighbor)")
    parser.add_option("model", "t'", float, 0.0, "Transfer integral (Second nearest)")
    parser.add_option("model", "nelec", float, 1.0, "Number of electrons per unit cell.")
    parser.add_option("model", "norb", str, "1",
                      "Number of orbitals at each inequivalent shell (integers separated by commas or spaces.)")
    parser.add_option("model", "ncor", int, 1, "Number of correlated shells in a unit cell (for lattice = wannier90).")
    parser.add_option("model", "corr_to_inequiv", str, "None",
                      "Mapping from correlated shells to equivalent shells (for lattice = wannier90)")
    parser.add_option("model", "bvec", str, "[(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)]", "Reciprocal lattice vectors in arbitrary unit.")
    parser.add_option("model", "nk", int, 8, "Number of *k* along each line")
    parser.add_option("model", "nk0", int, 0, "Number of *k* along b_0 (for lattice = wannier90, external)")
    parser.add_option("model", "nk1", int, 0, "Number of *k* along b_1 (for lattice = wannier90, external)")
    parser.add_option("model", "nk2", int, 0, "Number of *k* along b_2 (for lattice = wannier90, external)")
    parser.add_option("model", "spin_orbit", bool, False, "Whether the spin-orbit case.")
    parser.add_option("model", "interaction", str, "kanamori",
                      'Chosen from "slater_uj", "slater_f", "kanamori", "respack" (See below)')
    parser.add_option("model", "density_density", bool, False,
                      "If true, only the density-density part of the interaction is used (See below).")
    parser.add_option("model", "kanamori", str, "None",
                      "U (Diagonal Coulomb pot.), U\' (Off-diagonal Coulomb pot.) and J (Hund coupling) (See below).")
    parser.add_option("model", "slater_f", str, "None", "Angular momentum, Slater integrals F (See below).")
    parser.add_option("model", "slater_uj", str, "None", "Angular momentum, Slater integrals in U and J (See below).")
    parser.add_option("model", "local_potential_matrix", str, "None", "dict of {ish: 'filename'} to specify local potential matrix of ish-th shell")
    parser.add_option("model", "local_potential_factor", str, "1.0", "Prefactors to the local potential matrix (float or list with len=ncor)")

    # [system]
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature. This parameter is overridden, if T is given.")
    parser.add_option("system", "T", float, -1.0, "Temperature. If this parameter is given, beta is overridden by 1/T.")
    parser.add_option("system", "n_iw", int, 2048, "Number of Matsubara frequencies")
    parser.add_option("system", "fix_mu", bool, False, "Whether or not to fix chemical potential to a given value.")
    parser.add_option("system", "mu", float, 0.0, "Initial chemical potential.")
    parser.add_option("system", "prec_mu", float, 0.0001,
                      "Threshold for calculating chemical potential with the bisection method.")
    parser.add_option("system", "with_dc", bool, False, "Whether or not use double-counting correction (See below)")
    parser.add_option("system", "dc_type", str, 'HF_DFT', "Chosen from 'HF_DFT' (default), 'HF_imp', 'FLL'")

    # [impurity_solver]
    parser.add_option("impurity_solver", "name", str, 'null',
                      "Name of impurity solver. Available options are null, TRIQS/cthyb, TRIQS/hubbard-I, ALPS/cthyb, ALPS/cthyb-seg, pomerol.")
    parser.add_option("impurity_solver", "basis_rotation", str, 'None', "You can specify either 'Hloc', 'None', or the location of a file..")
    parser.allow_undefined_options("impurity_solver")

    # [control]
    parser.add_option("control", "max_step", int, 100, "Maximum steps of DMFT loops")
    parser.add_option("control", "sigma_mix", float, 0.5, "Mixing parameter for self-energy")
    parser.add_option("control", "restart", bool, False,
                      "Whether or not restart from a previous calculation stored in a HDF file.")
    parser.add_option("control", "initial_static_self_energy", str, "None", "dict of {ish: 'filename'} to specify initial value of the self-energy of ish-th shell. The file format is the same as local_potential_matrix.")
    parser.add_option("control", "initial_self_energy", str, "None", "Filename containing initial self-energy in the same format as sigma.dat generated by dcore_check.")
    parser.add_option("control", "time_reversal", bool, False, "If true, an average over spin components are taken.")
    parser.add_option("control", "symmetry_generators", str, 'None', "Generators for symmetrization of self-energy.")
    parser.add_option("control", "n_converge", int, 1, "The DMFT loop is terminated if the convergence criterion defined with converge_tol is satisfied n_converge times consecutively.")
    parser.add_option("control", "converge_tol", float, 0, "Tolerance in the convergence check. The chemical potential and the renormalization factor are examined.")

    # [tool]
    parser.add_option("tool", "nnode", int, 0, "[NOT USED] Number of node for the *k* path", OptionStatus.RETIRED)
    parser.add_option("tool", "nk_line", int, 8, "Number of *k* along each line")
    parser.add_option("tool", "knode", str, "[(G,0.0,0.0,0.0),(X,1.0,0.0,0.0)]",
                      "The name and the fractional coordinate of each k-node.")
    parser.add_option("tool", "omega_min", float, -1, "Minimum value of real frequency")
    parser.add_option("tool", "omega_max", float, 1, "Max value of real frequency")
    parser.add_option("tool", "Nomega", int, 100, "Number of real frequencies")
    parser.add_option("tool", "broadening", float, 0.0, "An additional Lorentzian broadening")
    parser.add_option("tool", "eta", float, 0.0, "Imaginary frequency shift for the Pade approximation")
    parser.add_option("tool", "omega_pade", float, 1E+20, "Cutoff frequency for the Pade approximation. Data in [-i omega_pade, i omega_pade] is used.")
    parser.add_option("tool", "n_pade_min", int, 20, "Minimum number of Matsubara frequencies used for Pade approximation.")
    parser.add_option("tool", "n_pade_max", int, -1, "Maximum number of Matsubara frequencies used for Pade approximation. If negative, this will be replaced with n_iw in [system] block.")

    parser.add_option("tool", "omega_check", float, 0, "Maximum frequency for dcore_check. If not specified, a fixed number of Matsubara points are taken.")
    parser.add_option("tool", "nk_mesh", int, 0, "Number of k points along each axis for computation of A(k,omega) on a 3D mesh")
    parser.add_option("tool", "nk0_mesh", int, 0, "Number of k points along b_0 for computation of A(k,omega) on a 3D mesh")
    parser.add_option("tool", "nk1_mesh", int, 0, "Number of k points along b_1 for computation of A(k,omega) on a 3D mesh")
    parser.add_option("tool", "nk2_mesh", int, 0, "Number of k points along b_2 for computation of A(k,omega) on a 3D mesh")

    # [bse]
    parser.add_option("bse", "num_wb", int, 1, "Number of bosonic frequencies (>0)")
    parser.add_option("bse", "num_wf", int, 10, "Number of fermionic frequencies (>0)")
    parser.add_option("bse", "h5_output_file", str, 'dmft_bse.h5', "Output HDF5 file for bse data")
    parser.add_option("bse", "skip_X0q_if_exists", bool, False, "Skip X_0(q) calc if file already exists")
    parser.add_option("bse", "skip_Xloc", bool, False, "Skip X_loc calc (for RPA)")
    parser.add_option("bse", "use_temp_file", bool, False, "Whether or not temporary file is used in computing X0_q. This option will reduce the memory footprints.")
    parser.add_option("bse", "X0q_qpoints_saved", str, 'quadrant', "Specifies for which q points X0q are saved in a HDF file. quadrant or path to a q_path.dat file.")


    return parser

def _cast_to_bool(p):
    if p is None:
        return False
    elif isinstance(p, str):
        return (p != 'None' and p != '')
    elif isinstance(p, bool):
        return p


def two_options_incompatible(params, option1, option2):
    if _cast_to_bool(params[option1[0]][option1[1]]) and _cast_to_bool(params[option2[0]][option2[1]]):
        raise RuntimeError("[{}][{}] and [{}][{}] cannot be True at the same time!".format(option1[0], option1[1], option2[0], option2[1]))


def parse_parameters(params):
    """
    Parse some parameters in a parameter

    :param params: dict (will be updated)
    :return:  None
    """

    if 'system' in params:
        if params['system']['T'] > 0:
            params['system']['beta'] = 1.0 / params['system']['T']
            params['system']['T'] = -1.0  # To make sure that only beta will be used in computations

    if 'control' in params:
        two_options_incompatible(params, ('control', 'restart'), ('control', 'initial_static_self_energy'))
        two_options_incompatible(params, ('control', 'initial_self_energy'), ('control', 'initial_static_self_energy'))

    if 'model' in params:
        ncor = params['model']['ncor']

        # Set [model][corr_to_inequiv] and [model][n_inequiv_shells]
        if params['model']['corr_to_inequiv'] == 'None':
            params['model']['corr_to_inequiv'] = numpy.arange(ncor)
            params['model']['n_inequiv_shells'] = ncor
        else:
            equiv_str_list = re.findall(r'[^\s,]+', params['model']['corr_to_inequiv'])
            corr_to_inequiv = numpy.array(list(map(int, equiv_str_list)))
            if len(corr_to_inequiv) != ncor:
                raise RuntimeError("Invalid number of elements in corr_to_inequiv!")
            params['model']['corr_to_inequiv'] = corr_to_inequiv
            params['model']['n_inequiv_shells'] = len(numpy.unique(corr_to_inequiv))
            if numpy.amin(corr_to_inequiv) < 0 or numpy.amin(corr_to_inequiv) > params['model']['n_inequiv_shells']:
                raise RuntimeError('Elements of corr_to_inequiv must be in the range of [0, n_inequiv_shells-1]!')

        # Set [model][norb_inequiv_sh]
        nsh = params['model']['n_inequiv_shells']
        params['model']['norb_inequiv_sh'] = numpy.array(list(map(int, re.findall(r'\d+', params['model']['norb']))))
        if len(params['model']['norb_inequiv_sh']) != nsh:
            raise RuntimeError("Wrong number of entries in norb!")

        # Set [model][norb_corr_sh]
        corr_to_inequiv = params['model']['corr_to_inequiv']
        params['model']['norb_corr_sh'] = numpy.array([params['model']['norb_inequiv_sh'][corr_to_inequiv[icrsh]] for icrsh in range(ncor)])

    if 'mpi' in params:
        # Expand enviroment variables
        params['mpi']['command'] = os.path.expandvars(params['mpi']['command'])

    if 'tool' in params:
        if params['tool']['n_pade_max'] < 0:
            params['tool']['n_pade_max'] = params['system']['n_iw']


def parse_knode(knode_string):
    """
    Parse knode

    Parameters
    ----------
    knode_string
        (label, k0, k1, k2) in the fractional coordinate

    Returns
    -------
    knode list of KNode

    """

    KNode = namedtuple('KNode', ('kvec', 'label'))

    knode_list = re.findall(r'\(\w+,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)', knode_string.replace(" ", ""))
    knode = []
    try:
        for _list in knode_list:
            _knode = [w for w in re.split(r'[)(,]', _list) if len(w) > 0]
            knode.append(KNode(label = _knode[0], kvec = numpy.array(list(map(float, _knode[1:4])))))
    except RuntimeError:
        raise RuntimeError("Error ! Format of knode is wrong.")
    return knode


def parse_bvec(bvec_string):
    """
    Parse bvec

    Parameters
    ----------
    bvec_string
        [(b0x, b0y, k0z),(b1x, b1y, k1z),(b2x, b2y, k2z)]

    Returns
    -------
    bvec numpy.ndarray shape=(3,3)
    Each "row" vector corresponds to a reciprocal vector.

    """

    #bvec_list = re.findall(r'\(\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)', p["model"]["bvec"])
    bvec_list = ast.literal_eval(bvec_string)
    if isinstance(bvec_list, list) and len(bvec_list) == 3:
        bvec = numpy.array(bvec_list, dtype=float)
        assert bvec.shape == (3,3)
    else:
        raise RuntimeError("Error ! Format of bvec is wrong.")
    return bvec
