

import numpy
from itertools import product
import sys
import h5py
import argparse

from irbasis_util.four_point_ph_view import FourPointPHView
from irbasis_util.tensor_regression import predict
from .tools import perform_fit, construct_prj, predict_xloc
from ..tools import mpi_split, float_to_complex_array, complex_to_float_array


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ph.py',
        description='tensor regression code.',
        epilog='end',
        usage='$ ',
        add_help=True)
    parser.add_argument('data_file', action='store', default=None, type=str, help="data_file")
    parser.add_argument('--bfreq', action='store', type=int, help='Bosonic frequency')
    parser.add_argument('--niter', default=50, type=int, help='Number of iterations')
    parser.add_argument('--D', action='store', type=int, help='Rank of decomposition')
    parser.add_argument('--Lambda', default=0, type=float, help='Lambda')
    parser.add_argument('--svcutoff', default=0, type=float, help='svcutoff')
    parser.add_argument('--num_wf', action='store', type=int, help='num of fermionic frequnencies for interpolation')
    parser.add_argument('--alpha', action='store', default=1e-5, type=float, help='alpha')
    parser.add_argument('--seed', action='store', default=100, type=int, help='seed')
    parser.add_argument('--rtol', action='store', default=1e-3, type=float, help='rtol')

    from mpi4py import MPI
    import os

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    is_master_node = comm.Get_rank() == 0

    args = parser.parse_args()
    if os.path.isfile(args.data_file) is False:
        print("data file does not exist.")
        sys.exit(-1)
    boson_freq = int(args.bfreq)
    D = args.D

    with h5py.File(args.data_file, 'r') as hf:
        beta = hf['/bse_sparse/args/beta'][()]
        freqs_PH_all = hf['/bse_sparse/freqs'][()]
        nfreqs_all = len(freqs_PH_all)

        # Count number of freqs for target bosonic frequency
        idx_tb = freqs_PH_all[:,0] == boson_freq
        n_freqs = numpy.sum(idx_tb)
        freqs_PH = freqs_PH_all[idx_tb,:]
        if n_freqs == 0:
            print("Data for boson_freq{} was not found.".format(boson_freq))
            sys.exit(-1)
        del freqs_PH_all

        # Prepare for distributing data over processes
        sizes, offsets = mpi_split(n_freqs, comm.size)
        n_freqs_local = sizes[rank]
        start, end = offsets[rank], offsets[rank]+sizes[rank]

        # Distribute data!
        xloc_keys = list(hf['/bse_sparse/0'].keys())
        xloc_local = numpy.zeros((len(xloc_keys), n_freqs_local), dtype=complex)
        for i, key in enumerate(xloc_keys):
            h5_key = '/bse_sparse/0/{}'.format(key)
            xloc_local[i, :] = float_to_complex_array(hf[h5_key][()])[idx_tb][start:end]
        num_o = len(xloc_keys)

    Lambda = args.Lambda
    wmax = Lambda / beta

    freqs_f = [tuple(freqs_PH[i, 1:]) for i in range(n_freqs)]
    freqs_f_local = numpy.array(freqs_f)[start:end, :]

    # Construct basis
    phb = FourPointPHView(boson_freq, Lambda, beta, args.svcutoff, True)
    prj = construct_prj(phb, freqs_f_local)

    # Fit
    xs = perform_fit(prj, xloc_local, D, args.niter, args.seed, args.alpha, args.rtol)

    xloc_local_fit = predict_xloc(prj, xs)
    xloc_abs_max = comm.allreduce(numpy.amax(numpy.abs(xloc_local)), op=MPI.MAX)
    xloc_abs_diff = comm.allreduce(numpy.amax(numpy.abs(xloc_local-xloc_local_fit)), op=MPI.MAX)

    # Save results into the HDF file
    if is_master_node:
        with h5py.File(args.data_file, 'a') as hf:
            prefix = '/bse_sparse/decomposed/0/{}/D{}'.format(boson_freq, D)
            if prefix in hf:
                del hf[prefix]
            for i, x in enumerate(xs):
                hf[prefix + '/x' + str(i)] = x
            for i, key in enumerate(xloc_keys):
                hf[prefix + '/keys/' + str(i)] = key

        print('')
        print('Xloc abs max = ', xloc_abs_max)
        print('Xloc max error in fit = ', xloc_abs_diff)

    # Interpolation in the 2D box of [-num_wf, num_wf-1]
    num_wf = args.num_wf
    n_freqs_box = (2*num_wf)**2
    sizes_box, offsets_box = mpi_split(n_freqs_box, comm.size)

    freqs_f_box = numpy.empty((n_freqs_box, 2), dtype=int)
    for idx, (i, j) in enumerate(product(list(range(2*num_wf)), repeat=2)):
        freqs_f_box[idx, :] = (i-num_wf, j-num_wf)
    freqs_f_box_local = freqs_f_box[offsets_box[rank]:offsets_box[rank]+sizes_box[rank], :]

    prj_box = construct_prj(phb, freqs_f_box_local)
    # xloc_box_local (frequency, orbital)
    xloc_box_local = predict(prj_box, xs)

    # Note: xloc_box_local (frequency, orbital)
    xloc_box = comm.gather(xloc_box_local, root=0)
    if is_master_node:
        # Join arrays along the frequency axis
        xloc_box = numpy.concatenate(xloc_box, axis=0)
        with h5py.File(args.data_file, 'a') as hf:
            prefix = '/bse_sparse/interpolated/0/wb{}/D{}'.format(boson_freq, D)
            if prefix in hf:
                del hf[prefix]
            for i, key in enumerate(xloc_keys):
                hf[prefix + '/' + key] = complex_to_float_array(numpy.array(xloc_box[:, i].reshape((2*num_wf, 2*num_wf))))
