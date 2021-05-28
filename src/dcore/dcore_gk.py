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
from warnings import warn

from h5 import HDFArchive

from .program_options import create_parser, parse_parameters
from .dmft_core import DMFTCoreSolver
from .sumkdft_workers.launcher import run_sumkdft
from .tools import complex_to_float_array

class _DMFTSolver(DMFTCoreSolver):
    def __init__(self, seedname, params, output_file='', output_group='dmft_out'):
        super().__init__(seedname, params, output_file, output_group, read_only=True, restart=True)

    def calc_gk(self, smpl_freqs):
        """
        Compute G(k, iw)
        """
        params = self._make_sumkdft_params()
        params['mu'] = self._chemical_potential

        # Number of positive fermionic frequencies
        n_iw = self._n_iw

        if numpy.abs(smpl_freqs).max() >= n_iw:
            warn("Some of sampling frequencies are not on the frequency mesh. "
                       "Data on these sampling points will be replaced by 1/iw."
                       "This may cause systematic errors if the frequency window size is small."
                       )
        lookup_smpl_freqs_win = numpy.where(numpy.abs(smpl_freqs) < n_iw)[0]
        params['smpl_freqs'] = smpl_freqs[lookup_smpl_freqs_win]
        r = run_sumkdft(
            'SumkDFTWorkerGk',
            os.path.abspath(self._seedname+'.h5'), './work/gk', self._mpirun_command, params)
        
        # Compute only for sampling frequencies within the window
        # Use 1/iw for those outside the window
        # nk, nfreqs, 2, norb, 2, norb
        gk_win = r['gk']
        nk, nfreqs_win, _, num_orb, _, _ = gk_win.shape
        num_so = 2*num_orb
        gk_win = gk_win.reshape((nk, nfreqs_win, 2*num_orb, 2*num_orb))
        gk = numpy.empty((nk, smpl_freqs.size, num_so, num_so), dtype=numpy.complex128)
        gk[:, lookup_smpl_freqs_win, :, :] = gk_win
        iw = 1J * (2*smpl_freqs + 1) * numpy.pi/self._beta

        # giw_tail: (n_smpl_freqs, num_so, num_so)
        giw_tail = numpy.einsum('w,ij->wij', 1/iw, numpy.identity(num_so))
        mask = numpy.abs(smpl_freqs) >= n_iw
        gk[:, mask, :, :] = giw_tail[None, mask, :, :]

        return gk

def dcore_gk(filename, np=1, prefix="./"):
    """
    Execute dcore_gk

    Parameters
    ----------
    filename : string
        Input-file name
    """
    print("\n############  Reading Input File  #################\n")
    print("  Input File Name : ", filename)

    # Construct a parser with default values
    pars = create_parser(['model', 'system', 'tool', 'mpi'])

    # Parse keywords and store
    pars.read(filename)
    p = pars.as_dict()
    parse_parameters(p)
    seedname = p["model"]["seedname"]
    p["mpi"]["num_processes"] = np

    ## make output directory
    dir = os.path.dirname(prefix)
    if not os.path.exists(dir):
        os.makedirs(dir)

    # Summary of input parameters
    print("\n  @ Parameter summary")
    print("\n    [model] block")
    for k, v in list(p["model"].items()):
        print("      {0} = {1}".format(k, v))
    print("\n    [tool] block")
    for k, v in list(p["tool"].items()):
        print("      {0} = {1}".format(k, v))

    #
    # Coompute G(k, iw)
    #
    print("\n#############   Computing G(k, iw) ####################\n")
    to_be_saved = {}
    if p["tool"]["gk_smpl_freqs"] == "sparse":
        import irbasis_x, irbasis
        b_xy = irbasis_x.twopoint.TruncatedBasis(irbasis.load('F', p['tool']['Lambda_IR']), cutoff=p['tool']['cutoff_IR'])
        smpl_freqs = b_xy.sampling_points_matsubara(b_xy.dim()-1)
        to_be_saved['Lambda_IR'] = p['tool']['Lambda_IR']
        to_be_saved['cutoff_IR'] = p['tool']['cutoff_IR']
    elif p["tool"]["gk_smpl_freqs"] == "dense":
        n_iw = p['system']['n_iw']
        smpl_freqs = numpy.arange(-n_iw, n_iw)
    else:
        with open(p["tool"]["gk_smpl_freqs"], 'r') as f:
            nfreqs = int(f.readline())
            smpl_freqs = []
            for i in range(nfreqs):
                elem = f.readline().split()
                if int(elem[0]) != i:
                    raise RuntimeError(f"""File format of {p["tool"]["gk_smpl_freqs"]} is wrong!""")
                smpl_freqs.append(int(elem[1]))
        smpl_freqs = numpy.array(smpl_freqs)


    solver = _DMFTSolver(seedname, p)
    gk = solver.calc_gk(smpl_freqs)
    nw = gk.shape[1]
    nf = gk.shape[-1]
    # (nk, w, f0, f1) => (w, f0, f1, nk)
    gk = numpy.moveaxis(gk, 0, -1)
    # (w, f0, f1, nk) => (w, spin0, orb0, spin1, orb1, nk0, nk1, nk2)
    gk = gk.reshape(
        (nw, 2, nf//2, 2, nf//2) +
        (p['model']['nk0'], p['model']['nk1'], p['model']['nk2'])
    )

    # Spin orbitals are ordered regardless of 
    # the value of "spin_orbit" option as follows:
    # for sh in correlated_shells:
    #     for spin in ['up', 'down']:
    with HDFArchive(os.path.abspath("{}./{}_gk.h5".format(prefix, seedname)), 'w') as h:
        h.create_group('gkw')
        for k, v in to_be_saved.items():
            h[k] = v
        h['beta'] = p['system']['beta']
        h['data'] = complex_to_float_array(gk)
        h['smpl_freqs'] = smpl_freqs

    print("\n#################  Done  #####################\n")


def run():
    import argparse
    from dcore.option_tables import generate_all_description
    from dcore.version import version

    parser = argparse.ArgumentParser(
        prog='dcore_gk.py',
        description='post-processing script for dcore.',
        usage='$ dcore_gk input --np 4',
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
    parser.add_argument('--np', help='Number of MPI processes', required=True)
    parser.add_argument('--version', action='version', version='DCore {}'.format(version))
    parser.add_argument('--prefix',
                        action='store',
                        default='./',
                        type=str,
                        help='prefix for output files (default: ./)'
                        )

    args = parser.parse_args()
    if os.path.isfile(args.path_input_file) is False:
        print("Input file does not exist.")
        sys.exit(-1)
    dcore_gk(args.path_input_file, int(args.np), args.prefix)
