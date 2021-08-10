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

from h5 import HDFArchive
from ..tools import launch_mpi_subprocesses, raise_if_mpi_imported, make_empty_dir

def run_sumkdft(runner_cls, model_file, work_dir, mpirun_command, params):
    """
    Runs SumKDFT by launching MPI processes.

    :param runner_cls: str
       Name of subclass of SumkDFTWorker
    :param model_file: str
        HDF5 file
    :param work_dir: str
        Working directory
    :param mpirun_command: str
        Command for executing mpi programs
    :param params: dict
        Parameters for SumkDFT
    :return: dict
        results

    params contains the following parameters.
        mu            : float, chemical potential. If mu is not given, mu will be adjusted (optional).
        prec_mu       : float, precision of adjustment of chemical potential (optional)
        broadening    : float, broadening parameter for DOS (must be set when calc_mode = dos, spaghettis)
        mesh          : (float, float, int) real-frequency mesh (optional)

    """

    raise_if_mpi_imported()

    # Prepare input files
    make_empty_dir(work_dir)

    cwd_org = os.getcwd()
    os.chdir(work_dir)

    if os.path.exists('./input.h5'):
        os.remove('./input.h5')
    if os.path.exists('./output.h5'):
        os.remove('./output.h5')

    with HDFArchive('./input.h5', 'w') as h:
        h['params'] = params

    commands = [sys.executable, "-m", "dcore.sumkdft_workers.mpi_main"]
    commands.append(runner_cls)
    commands.append(model_file)
    commands.append(os.path.abspath('./input.h5'))
    commands.append(os.path.abspath('./output.h5'))

    with open('./output', 'w') as output_file:
        launch_mpi_subprocesses(mpirun_command, commands, output_file)

    with open('./output', 'r') as output_file:
        for line in output_file:
            print(line, end='')

    results = {}
    if os.path.exists(os.path.abspath('./output.h5')):
        with HDFArchive(os.path.abspath('./output.h5'), 'r') as h:
            for k in list(h.keys()):
                results[k] = h[k]

    os.chdir(cwd_org)

    return results