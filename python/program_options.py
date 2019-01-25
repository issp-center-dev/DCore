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


def create_parser():
    """
    Create a parser for all program options of DCore
    """
    parser = TypedParser()

    # [model]
    parser.add_option("model", "t", float, 1.0, "Transfer integral (Nearest neighbor)")
    parser.add_option("model", "t'", float, 0.0, "Transfer integral (Second nearest)")
    parser.add_option("model", "ncor", int, 1, "Number of correlated shells in a unit cell (Only wannier90).")
    parser.add_option("model", "lattice", str, "chain",
                      'Chosen from "chain", "square", "cubic", "bethe", and "wannier90"')
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
                      "If true, only the density-density part of the interaction is used (See below).")
    parser.add_option("model", "kanamori", str, "None",
                      "U (Diagonal Coulomb pot.), U\' (Off-diagonal Coulomb pot.) and J (Hund coupling) (See below).")
    parser.add_option("model", "slater_f", str, "None", "Angular momentum, Slater integrals F (See below).")
    parser.add_option("model", "slater_uj", str, "None", "Angular momentum, Slater integrals in U and J (See below).")
    parser.add_option("model", "non_colinear", bool, False,
                      "Set True for the case that non-colinear DMFT from the COLINEAR DFT calculation.")

    # [system]
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature.")
    parser.add_option("system", "n_iw", int, 2048, "Number of Matsubara frequencies")
    parser.add_option("system", "n_tau", int, 10000, "Number of imaginary-time points")
    parser.add_option("system", "fix_mu", bool, False, "Whether or not to fix chemical potential to a given value.")
    parser.add_option("system", "mu", float, 0.0, "Initial chemical potential.")
    parser.add_option("system", "nk", int, 8, "Number of *k* along each line")
    parser.add_option("system", "nk0", int, 0, "Number of *k* along b_0 (only for wannier90)")
    parser.add_option("system", "nk1", int, 0, "Number of *k* along b_1 (only for wannier90)")
    parser.add_option("system", "nk2", int, 0, "Number of *k* along b_2 (only for wannier90)")
    parser.add_option("system", "prec_mu", float, 0.0001,
                      "Threshold for calculating chemical potential with the bisection method.")
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature.")
    parser.add_option("system", "with_dc", bool, False, "Whether or not use double counting correction (See below)")
    parser.add_option("system", "perform_tail_fit", bool, False, "Whether or not perform the tail-fit.")
    parser.add_option("system", "fit_max_moment", int, 2, "Highest moment to fit in the tail of Sigma_iw.")
    parser.add_option("system", "fit_min_w", float, 5.0, "Matsubara frequency from which tail fitting should start.")
    parser.add_option("system", "fit_max_w", float, 10.0, "Matsubara frequency at which tail fitting should end.")
    parser.add_option("system", "n_l", int, 0,
                      "The number of the Legendre polynomial for QMC. If not, the solver's default value is used.")

    # [impurity_solver]
    parser.add_option("impurity_solver", "name", str, 'TRIQS/hubbard-I',
                      "Name of impurity solver. Available options are TRIQS/cthyb, TRIQS/hubbard-I, ALPS/cthyb.")
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

    return parser
