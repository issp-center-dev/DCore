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
import os
import sys
import copy
import irbasis_x
from irbasis_x import bse_dmft

from .dmft_core import DMFTCoreSolver
from .program_options import create_parser, parse_parameters
from .tools import *
import h5py
from .irbasis_util import construct_basis
from . import impurity_solvers

def calc_Floc_impurity_model(solver_name, solver_params, mpirun_command, basis_rot,
                              Umat, gf_struct, beta, n_iw,
                              Sigma_iw, Gloc_iw, ish, wsample_ph):
    """
    Calculate Floc of an impurity model
    """
    wsample_ph = irbasis_x.freq.check_ph_convention(*wsample_ph)

    Solver = impurity_solvers.solver_classes[solver_name]

    ##if not Solver.is_Floc_computable():
        #raise RuntimeError(f"Floc is not computable with {solver_name}!")

    raise_if_mpi_imported()

    sol = Solver(beta, gf_struct, Umat, n_iw)

    G0_iw = dyson(Sigma_iw=Sigma_iw, G_iw=Gloc_iw)
    sol.set_G0_iw(G0_iw)

    # Compute rotation matrix to the diagonal basis if supported
    rot = impurity_solvers.compute_basis_rot(basis_rot, sol)
    s_params = copy.deepcopy(solver_params)
    s_params['random_seed_offset'] = 1000 * ish

    work_dir_org = os.getcwd()
    work_dir = 'work/imp_shell' + str(ish) + '_bse'
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)

    # Solve the model
    if Solver.is_Floc_computable():
        Floc = sol.calc_Floc_ph_sparse(rot, mpirun_command, wsample_ph, s_params)
        G2loc = None
    else:
        Floc = None
        G2loc = sol.calc_G2loc_ph_sparse(rot, mpirun_command, wsample_ph, s_params)

    os.chdir(work_dir_org)

    return Floc, G2loc


class DMFTBSESolver(DMFTCoreSolver):
    def __init__(self, seedname, params, output_file='', output_group='dmft_out'):
        super().__init__(seedname, params, output_file, output_group, read_only=True, restart=True)

    def compute_G2loc(self):
        """
        Compute local two-particle gf
        """
        p = self._params

        solver_name = p['impurity_solver']['name']
        Lambda_IR = p['tool']['Lambda_IR']
        cutoff_IR = p['tool']['cutoff_IR_2P']

        # Sampling bosonic frequencies (i.e, -1, 0, 1, 2, ...)
        wb_str = p['vertex']['wb_ph']
        if wb_str == '':
            raise RuntimeError("Please set wb_ph!")
        elif wb_str == 'sparse':
            basis_b = construct_basis('B', p['system']['beta'], Lambda_IR, cutoff_IR)
            wb_sample = basis_b.wsample//2
            wb_sample = wb_sample[wb_sample >= 0]
        else:
            # Try to interpret wb_str as integers
            str_ = wb_str.replace(',', ' ').split()
            wb_sample = numpy.array(list(map(int, str_)))

        # Compute local Green's function
        Gloc_iw_sh, _ = self.calc_Gloc()

        # Construct a DMFT BSE solver
        #   Note: In irbasis_x, fermionic/bosonic frequences are denoted by odd/even int numbers.
        basis = irbasis_x.FourPointBasis(Lambda_IR, self._beta, cutoff_IR)
        wsample_ph = [], [], []
        for wb_ in wb_sample:
            solver = bse_dmft.SparseSolver(basis, 2*wb_)
            for i in range(3):
                wsample_ph[i].append(solver.wsample_Floc[i])
        wsample_ph = tuple(map(numpy.hstack, wsample_ph))

        with h5py.File(self._seedname + '_vertex.h5', 'w') as h:
            h['wb_sample'] = wb_sample
            h['Lambda_IR'] = Lambda_IR
            h['cutoff_IR'] = cutoff_IR
            h['n_inequiv_sh'] = self._n_inequiv_shells
            h['corr_to_inequiv'] = numpy.array(self.corr_to_inequiv)
            if self.use_spin_orbit:
                h['nso_sh'] = numpy.array(self._dim_sh)
            else:
                h['nso_sh'] = 2 * numpy.array(self._dim_sh)
            h['wsample_ph'] = wsample_ph
            for ish in range(self._n_inequiv_shells):
                print("\nSolving impurity model for inequivalent shell " + str(ish) + " ...")
                sys.stdout.flush()
                Floc, G2loc = calc_Floc_impurity_model(
                    solver_name, self._solver_params, self._mpirun_command,
                    self._params["impurity_solver"]["basis_rotation"],
                    self._Umat[ish], self._gf_struct[ish],
                    self._beta, self._n_iw,
                    self._sh_quant[ish].Sigma_iw, Gloc_iw_sh[ish],
                    ish, wsample_ph)
                if Floc is not None:
                    Floc = numpy.moveaxis(Floc, -1, 0) # (nfreq, nf, nf, nf, nf)
                    nw, nf = Floc.shape[0], Floc.shape[1]
                    Floc = Floc.reshape((nw,) + (2, nf//2)* 4)
                    h[f'Floc/sh{ish}'] = complex_to_float_array(Floc)
                else:
                    assert G2loc is not None
                    G2loc = numpy.moveaxis(G2loc, -1, 0) # (nfreq, nf, nf, nf, nf)
                    nw, nf = G2loc.shape[0], G2loc.shape[1]
                    G2loc = G2loc.reshape((nw,) + (2, nf//2)* 4)
                    h[f'G2loc/sh{ish}'] = complex_to_float_array(G2loc)


def dcore_vertex(filename, np=1):
    """
    Main routine for the BSE post-processing tool

    Parameters
    ----------
    filename : string
        Input-file name
    """

    print("\n############  Reading Input File  #################\n")
    print("  Input File Name : ", filename)

    #
    # Construct a parser with default values
    #
    pars = create_parser(['model', 'system', 'impurity_solver', 'mpi', 'tool', 'vertex'])

    #
    # Parse keywords and store
    #
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    seedname = p["model"]["seedname"]
    p["mpi"]["num_processes"] = np

    # Load DMFT data
    solver = DMFTBSESolver(seedname, p, output_file=seedname + '.out.h5')
    if solver.iteration_number == 0:
        raise RuntimeError("Number of iterations is zero!")
    print("Number of iterations :", solver.iteration_number)

    # Compute and save vertex functions
    solver.compute_G2loc()

    # Finish
    print("\n#################  Done  #####################\n")


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version

    parser = argparse.ArgumentParser(
        prog='dcore_sparse_vertex.py',
        description='Compute vertex functions',
        usage='$ dcore_vertex input --np 4',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=generate_all_description(),
        add_help=True)
    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name."
                        )
    parser.add_argument('--np', help='Number of MPI processes', required=True)
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))
    #parser.add_argument('--mode', help='Calculation mode (calc_X0q, calc_G2loc, solve_BSE)', type=str, required=True)

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)

    dcore_vertex(args.path_input_file, int(args.np))