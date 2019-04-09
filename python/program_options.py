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

from typed_parser import *

import numpy
import re

def create_parser():
    """
    Create a parser for all program options of DCore
    """
    parser = TypedParser()

    # [mpi]
    parser.add_option("mpi", "command", str, "mpirun -np #", "Command for executing a MPI job. # will be relaced by the number of processes.")

    # [model]
    parser.add_option("model", "t", float, 1.0, "Transfer integral (Nearest neighbor)")
    parser.add_option("model", "t'", float, 0.0, "Transfer integral (Second nearest)")
    parser.add_option("model", "ncor", int, 1, "Number of correlated shells in a unit cell (Only wannier90).")
    parser.add_option("model", "lattice", str, "chain",
                      'Chosen from "chain", "square", "cubic", "bethe", "wannier90", and "external"')
    parser.add_option("model", "nelec", float, 1.0, "Number of electrons per unit cell.")
    parser.add_option("model", "seedname", str, "dcore", "Name of the system. The model HDF5 file will be seedname.h5.")
    parser.add_option("model", "norb", str, "1",
                      "Number of orbitals at each correlated shell (*ncor* integers separated by commas or spaces.)")
    parser.add_option("model", "equiv", str, "None",
                      "Equivalence of each correlated shell. Please, be careful to use it (See below).")
    parser.add_option("model", "bvec", str, "[(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)]", "Reciprocal lattice vectors in arbitrary unit.")
    parser.add_option("model", "spin_orbit", bool, False, "Whether the spin-orbit case (See :ref:`pbtutorial`).")
    parser.add_option("model", "time_reversal", bool, False, "If true, an average over spin components are taken.")
    parser.add_option("model", "interaction", str, "kanamori",
                      'Chosen from "slater_uj", "slater_f", "kanamori", "respack" (See below)')
    parser.add_option("model", "density_density", bool, False,
                      "If true, only the density-density part of the interaction is used (See below).",
                      OptionStatus.RETIRED)
    parser.add_option("model", "kanamori", str, "None",
                      "U (Diagonal Coulomb pot.), U\' (Off-diagonal Coulomb pot.) and J (Hund coupling) (See below).")
    parser.add_option("model", "slater_f", str, "None", "Angular momentum, Slater integrals F (See below).")
    parser.add_option("model", "slater_uj", str, "None", "Angular momentum, Slater integrals in U and J (See below).")
    parser.add_option("model", "non_colinear", bool, False,
                      "Set True for the case that non-colinear DMFT from the COLINEAR DFT calculation.")
    parser.add_option("model", "local_potential_matrix", str, "None", "dict of {ish: 'filename'} to specify local potential matrix of ish-th shell")
    parser.add_option("model", "local_potential_factor", str, "1.0", "Prefactors to the local potential matrix (float or list with len=ncor)")

    # [system]
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature.")
    parser.add_option("system", "n_iw", int, 2048, "Number of Matsubara frequencies")
    parser.add_option("system", "n_tau", int, 10000, "Number of imaginary-time points")
    parser.add_option("system", "fix_mu", bool, False, "Whether or not to fix chemical potential to a given value.")
    parser.add_option("system", "mu", float, 0.0, "Initial chemical potential.")
    parser.add_option("system", "nk", int, 8, "Number of *k* along each line")
    parser.add_option("system", "nk0", int, 0, "Number of *k* along b_0 (for lattice = wannier90, external)")
    parser.add_option("system", "nk1", int, 0, "Number of *k* along b_1 (for lattice = wannier90, external)")
    parser.add_option("system", "nk2", int, 0, "Number of *k* along b_2 (for lattice = wannier90, external)")
    parser.add_option("system", "prec_mu", float, 0.0001,
                      "Threshold for calculating chemical potential with the bisection method.")
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature.")
    parser.add_option("system", "with_dc", bool, False, "Whether or not use double counting correction (See below)")
    parser.add_option("system", "perform_tail_fit", bool, False, "Whether or not perform the tail-fit.")
    parser.add_option("system", "fit_max_moment", int, 2, "Highest moment to fit in the tail of Sigma_iw.")
    parser.add_option("system", "fit_min_w", float, 5.0, "Matsubara frequency from which tail fitting should start.")
    parser.add_option("system", "fit_max_w", float, 10.0, "Matsubara frequency at which tail fitting should end.")
    parser.add_option("system", "n_l", int, 0,
                      "The number of the Legendre polynomial for QMC. If not, the solver's default value is used.",
                      OptionStatus.RETIRED)
    parser.add_option("system", "initial_self_energy", str, "None", "dict of {ish: 'filename'} to specify initial value of the self-energy of ish-th shell. The file format is the same as local_potential_matrix.")

    # [impurity_solver]
    parser.add_option("impurity_solver", "name", str, 'null',
                      "Name of impurity solver. Available options are null, TRIQS/cthyb, TRIQS/hubbard-I, ALPS/cthyb, ALPS/cthyb-seg.")
    parser.add_option("impurity_solver", "basis_rotation", bool, True,
                      "If True, impurity models are solved in the basis that diagonalizes the local non-interacting Hamiltonian.")
    parser.allow_undefined_options("impurity_solver")

    # [control]
    parser.add_option("control", "max_step", int, 100, "Maximum steps of DMFT loops")
    parser.add_option("control", "sigma_mix", float, 0.5, "Mixing parameter for self-energy")
    parser.add_option("control", "restart", bool, False,
                      "Whether or not restart from a previous calculation stored in a HDF file.")

    # [tool]
    parser.add_option("tool", "nnode", int, 0, "[NOT USED] Number of node for the *k* path", OptionStatus.RETIRED)
    parser.add_option("tool", "nk_line", int, 8, "Number of *k* along each line")
    parser.add_option("tool", "knode", str, "[(G,0.0,0.0,0.0),(X,1.0,0.0,0.0)]",
                      "The name and the fractional coordinate of each k-node.")
    parser.add_option("tool", "omega_min", float, -1, "Minimum value of real frequency")
    parser.add_option("tool", "omega_max", float, 1, "Max value of real frequency")
    parser.add_option("tool", "Nomega", int, 100, "Number of real frequencies")
    parser.add_option("tool", "broadening", float, 0.1, "An additional Lorentzian broadening")
    parser.add_option("tool", "eta", float, 0.0, "Imaginary frequency shift for the Pade approximation")
    parser.add_option("tool", "omega_pade", float, 5.0, "Cutoff frequencies for the Pade approximation")
    parser.add_option("tool", "omega_check", float, 5.0, "Maximum frequency for dcore_check.")

    # [bse]
    parser.add_option("bse", "num_wb", int, 0, "Number of bosonic frequencies (>=0)")
    parser.add_option("bse", "num_wf", int, 10, "Number of fermionic frequencies (>0)")
    parser.add_option("bse", "h5_output_file", str, 'dmft_bse.h5', "Output HDF5 file for bse data")
    parser.add_option("bse", "skip_X0q_if_exists", bool, False, "Skip X_0(q) calc if file already exists")
    parser.add_option("bse", "skip_Xloc", bool, False, "Skip X_loc calc (for RPA)")
    parser.add_option("bse", "use_temp_file", bool, False, "Whether or not temporary file is used in computing X0_q. This option will reduce the memory footprints.")
    parser.add_option("bse", "X0q_qpoints_saved", str, 'quadrant', "Specifies for which q points X0q are saved in a HDF file. quadrant or path to a q_path.dat file.")

    return parser


def parse_parameters(params):
    """
    Parse some parameters in a parameter

    :param params: dict (will be updated)
    :return:  None
    """

    params['model']['norb_corr_sh'] = numpy.array(map(int, re.findall(r'\d+', params['model']['norb'])))

    ncor = params['model']['ncor']
    if params['model']['equiv'] == 'None':
        params['model']['equiv_sh'] = numpy.arange(ncor)
    else:
        equiv_list = re.findall(r'[^\s,]+', params['model']['equiv'])
        equiv_str_list = []
        equiv_index = 0
        equiv = numpy.zeros(ncor, dtype=int)
        for icor in range(ncor):
            if equiv_list[icor] in equiv_str_list:
                # Defined before
                equiv[icor] = equiv_str_list.index(equiv_list[icor])
            else:
                # New one
                equiv_str_list.append(equiv_list[icor])
                equiv[icor] = equiv_index
                equiv_index += 1
        params['model']['equiv_sh'] = equiv
