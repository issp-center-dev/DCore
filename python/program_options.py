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
    parser.add_option("model", "U", float, 0.0, "On-site Coulomb potential")
    parser.add_option("model", "J", float, 0.0, "On-site Hund potential")
    parser.add_option("model", "orbital_model", str, "single", 'Chosen from "single", "eg", "t2g", "full-d"')
    parser.add_option("model", "ncor", int, 1, "Number of correlation shell (Only wannier90).")
    parser.add_option("model", "lattice", str, "chain", 'Chosen from "chain", "square", "cubic", "bethe", and "wannier90"')
    parser.add_option("model", "nelec", float, 1.0, "Number of electrons per unit cell.")
    parser.add_option("model", "seedname", str, "dcore", "Name of the system. The model HDF5 file will be seedname.h5.")
    parser.add_option("model", "cshell", str, "[]", "Anguler momentum, and the number of orbitals of each correlation shell (Only wannier90). If not explicitly specified, the default value will be  [(0,1),...].")
    parser.add_option("model", "bvec", str, "[(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)]", "Reciprocal lattice vectors")

    # [system]
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature.")
    parser.add_option("system", "n_iw", int, 2048, "Number of Matsubara frequencies")
    parser.add_option("system", "n_tau", int, 10000, "Number of imaginary-time points")
    parser.add_option("system", "dc_type", int, -1, "Type of double-counting correction. Set -1 to disable DC correction. Other availale options are according to TRIQS DFTTools.")
    parser.add_option("system", "fix_mu", bool, False, "Whether or not to fix chemical potential to a given value.")
    parser.add_option("system", "mu", float, 1E+8, "Chemical potential used when fix_mu = True")
    parser.add_option("system", "nk", int, 8, "Number of *k* along each line")
    parser.add_option("system", "nk0", int, 0, "Number of *k* along b_0 (only for wannier90)")
    parser.add_option("system", "nk1", int, 0, "Number of *k* along b_1 (only for wannier90)")
    parser.add_option("system", "nk2", int, 0, "Number of *k* along b_2 (only for wannier90)")
    parser.add_option("system", "prec_mu", float, 0.0001, "Threshold for calculating chemical potential with the bisection method.")
    parser.add_option("system", "beta", float, 1.0, "Inverse temperature.")

    # [impurity_solver] 
    parser.add_option("impurity_solver", "name", str, 'TRIQS/hubbard-I', "Name of impurity solver. Available options are TRIQS/cthyb, TRIQS/hubbard-I, ALPS/cthyb.")
    parser.allow_undefined_options("impurity_solver")

    # [control] 
    parser.add_option("control", "max_step", int, 100, "Maximum steps of DMFT loops")
    parser.add_option("control", "sigma_mix", float, 0.5, "Mixing parameter for self-energy")
    parser.add_option("control", "delta_mix", float, 0.5, "Mixing parameter for hybridization function")
    parser.add_option("control", "restart", bool, False, "Whether or not restart from a previous calculation stored in a HDF file.")

    # [tool] 
    parser.add_option("tool", "nnode", int, 2, "Number of node for the *k* path")
    parser.add_option("tool", "nk_line", int, 8, "Number of *k* along each line")
    parser.add_option("tool", "knode", str, "[(G,0.0,0.0,0.0),(X,1.0,0.0,0.0)]", "The name and the fractional coordinate of each k-node.")
    parser.add_option("tool", "omega_min", float, -1, "Minimum value of real frequency")
    parser.add_option("tool", "omega_max", float, 1, "Max value of real frequency")
    parser.add_option("tool", "Nomega", int, 100, "Number of real frequencies")
    parser.add_option("tool", "broadening", float, 0.1, "An additional Lorentzian broadening")
    parser.add_option("tool", "eta", float, 0.0, "Imaginary frequency shift for the Pade approximation")
    parser.add_option("tool", "n_pade", int, 100, "Number of imaginary frequencies for the Pade approximation")

    return parser
