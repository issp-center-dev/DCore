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

import argparse
import os.path
import numpy
from dcore._dispatcher import MeshReFreq
from dcore.version import version, print_header
from dcore.program_options import create_parser, parse_parameters
from dcore.anacont import pade, spm

def dcore_anacont(inifile):
    parser = create_parser(["system", "model", "post", "post.anacont"])
    parser.read(inifile)
    params = parser.as_dict()
    parse_parameters(params)

    beta = params["system"]["beta"]

    omega_min = params["post"]["omega_min"]
    omega_max = params["post"]["omega_max"]
    if omega_min >= omega_max:
        # ToDo: stop the program properly
        assert omega_min < omega_max

    Nomega = params["post"]["Nomega"]
    mesh_w = MeshReFreq(omega_min, omega_max, Nomega)

    seedname = params["model"]["seedname"]
    file_sigma_iw = seedname + "_sigma_iw.npz"
    # file_sigma_iw = params["post"]["file_sigma_iw"]
    if not os.path.exists(file_sigma_iw):
        assert False, "File not found: " + file_sigma_iw
    sigma_iw_npz = numpy.load(file_sigma_iw)

    dir_work = params["post"]["dir_work"]
    if not os.path.exists(dir_work):
        os.makedirs(dir_work)

    solver = params["post.anacont"]["solver"]
    if solver == "pade":
        params_ac = pade.parameters_from_ini(inifile)
        data_w = pade.anacont(sigma_iw_npz, beta, mesh_w, params_ac)
    elif solver == "spm":
        params_ac = spm.parameters_from_ini(inifile)
        data_w = spm.anacont(sigma_iw_npz, beta, mesh_w, params_ac)
    else:
        assert False, "Unknown solver: " + solver

    file_sigma_w = os.path.join(dir_work, "sigma_w.npz")
    print("Writing to", file_sigma_w + "...")
    numpy.savez(file_sigma_w, **data_w)

def run():

    print_header()

    parser = argparse.ArgumentParser(
        prog='dcore_anacont.py',
        description='DCore post script -- analytic continuation.',
        usage='$ dcore_anacont input.ini',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        #epilog=generate_all_description()
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

    dcore_anacont(args.path_input_files)