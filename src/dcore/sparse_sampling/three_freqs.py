

import numpy
from itertools import product
import sys
import h5py
import argparse

from irbasis_util.four_point import FourPoint, from_PH_convention
from .tools import fit, construct_prj_three_freqs, predict_xloc
from ..tools import mpi_split, float_to_complex_array, complex_to_float_array


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='three_freqs.py',
        description='tensor regression code.',
        epilog='end',
        usage='$ ',
        add_help=True)
    parser.add_argument('data_file', action='store', default=None, type=str, help="data_file")
    parser.add_argument('--niter', default=50, type=int, help='Number of iterations')
    parser.add_argument('--D', action='store', type=int, help='Rank of decomposition')
    parser.add_argument('--Lambda', default=0, type=float, help='Lambda')
    parser.add_argument('--svcutoff', default=0, type=float, help='svcutoff')
    parser.add_argument('--num_wf', action='store', type=int, help='num of fermionic frequnencies for interpolation')
    parser.add_argument('--num_wb', action='store', type=int, help='num of bosonic frequnencies for interpolation')

    from mpi4py import MPI
    import os

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    is_master_node = comm.Get_rank() == 0

    args = parser.parse_args()
    if os.path.isfile(args.data_file) is False:
        print("data file is not exist.")
        sys.exit(-1)
    D = args.D

    with h5py.File(args.data_file, 'r') as hf:
        beta = hf['/bse_sparse/args/beta'][()]
        freqs_PH = hf['/bse_sparse/freqs'][()]
        n_freqs = len(freqs_PH)

        # Prepare for distributing data over processes
        sizes, offsets = mpi_split(n_freqs, comm.size)
        n_freqs_local = sizes[rank]
        start, end = offsets[rank], offsets[rank]+sizes[rank]

        # Distribute data!
        xloc_keys = list(hf['/bse_sparse/0'].keys())
        xloc_local = numpy.zeros((len(xloc_keys), n_freqs_local), dtype=complex)
        for i, key in enumerate(xloc_keys):
            h5_key = '/bse_sparse/0/{}'.format(key)
            xloc_local[i, :] = float_to_complex_array(hf[h5_key][()])[start:end]
        num_o = len(xloc_keys)

    Lambda = args.Lambda
    wmax = Lambda / beta

    freqs_PH = [tuple(freqs_PH[i, 1:]) for i in range(n_freqs)]
    freqs_PH_local = numpy.array(freqs_PH)[start:end, :]

    # Construct basis
    basis = FourPoint(Lambda, beta, args.svcutoff, True)
    prj = construct_prj_three_freqs(basis, freqs_PH_local)

    # Fit
    xs, mse = fit(prj, xloc_local, D, args.niter)

    # FIXME: ORDER OF FREQ AND ORBITAL
    xloc_local_fit = predict_xloc(prj, xs)
    xloc_abs_max = comm.allreduce(numpy.amax(numpy.abs(xloc_local)), op=MPI.MAX)
    xloc_abs_diff = comm.allreduce(numpy.amax(numpy.abs(xloc_local-xloc_local_fit)), op=MPI.MAX)

    # Save results into the HDF file
    if is_master_node:
        with h5py.File(args.data_file, 'a') as hf:
            prefix = '/bse_sparse/decomposed_three_freqs/0/D{}'.format(D)
            if prefix in hf:
                del hf[prefix]
            for i, x in enumerate(xs):
                hf[prefix + '/x' + str(i)] = x
            for i, key in enumerate(xloc_keys):
                hf[prefix + '/keys/' + str(i)] = key

        print('')
        print('Xloc abs max = ', xloc_abs_max)
        print('Xloc mean squared error in fit = ', mse)
        print('Xloc max error in fit = ', xloc_abs_diff)

    # Interpolation in the 3D box of [-num_wf:num_wf, -num_wf:num_wf, 0:num_wb]
    num_wf = args.num_wf
    num_wb = args.num_wb
    n_freqs_box = (2*num_wf)**2
    sizes_box, offsets_box = mpi_split(n_freqs_box, comm.size)
    for boson_freq in range(num_wb):
        # List of frequencies in the notation of four fermionic frequencies
        freqs_box = numpy.empty((n_freqs_box, 4), dtype=int)
        for idx, (i, j) in enumerate(product(list(range(2*num_wf)), repeat=2)):
            freqs_box[idx, :] = from_PH_convention((i-num_wf, j-num_wf, boson_freq))
        freqs_box_local = freqs_box[offsets_box[rank]:offsets_box[rank]+sizes_box[rank], :]

        prj_box = construct_prj_three_freqs(basis, freqs_box_local)
        xloc_box_local = predict_xloc(prj_box, xs)

        # Note: xloc_box_local (orbital, frequency)
        xloc_box = comm.gather(xloc_box_local, root=0)
        if is_master_node:
            # Join arrays along the frequency axis
            xloc_box = numpy.concatenate(xloc_box, axis=1).reshape((-1, num_wf, num_wf))
            with h5py.File(args.data_file, 'a') as hf:
                prefix = '/bse_sparse/interpolated/0/wb{}/D{}'.format(boson_freq, D)
                if prefix in hf:
                    del hf[prefix]
                for i, key in enumerate(xloc_keys):
                    hf[prefix + '/' + key] = complex_to_float_array(numpy.array(xloc_box[i, :, :]))
