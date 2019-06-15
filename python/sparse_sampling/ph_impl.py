from __future__ import print_function

import numpy
import copy

from irbasis_util.four_point_ph_view import FourPointPHView
from irbasis_util.tensor_regression import enable_MPI, optimize_als, OvercompleteGFModel, predict

from mpi4py import MPI 

enable_MPI() #for tensor_regression

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
is_master_node = comm.Get_rank() == 0

def construct_prj(phb, freqs_f_list):
    assert isinstance(phb, FourPointPHView)

    assert isinstance(freqs_f_list, numpy.ndarray)
    assert freqs_f_list.dtype == int
    assert freqs_f_list.ndim == 2
    assert freqs_f_list.shape[1] == 2

    prj = phb.projector_to_matsubara_vec(freqs_f_list, decomposed_form=True)
    for i in range(len(prj)):
        prj[i] = prj[i].reshape((-1, 12, phb.Nl))
    return prj

def fit(projectors, xloc_local, D, niter):
    num_o, num_freqs = xloc_local.shape
    num_rep = 12
    linear_dim = projectors[0].shape[-1] # dim of IR basis

    # (orbital, freq) => (freq, orbital)
    y = xloc_local.transpose((1,0))
    alpha_init = 0
    model = OvercompleteGFModel(num_freqs, num_rep, 2, num_o, linear_dim, projectors, y, alpha_init, D)
    info = optimize_als(model, niter, rtol=1e-20, optimize_alpha=1e-5, verbose=1, print_interval=1)

    return copy.deepcopy(model.x_tensors()), info['rmses'][-1]**2

def predict_xloc(prj, x_tensors):
    return predict(prj, x_tensors).transpose()


def interpolate(phb, freqs_f_list, x_tensors):
    assert isinstance(phb, FourPointPHView)
    assert isinstance(freqs_f_list, list)
    assert isinstance(x_tensors, list)

    n_freqs = len(freqs_f_list)
    prj = phb.projector_to_matsubara_vec(freqs_f_list, decomposed_form=True)
    for p in prj:
        print(p.shape)
    #for i in range(2):
    #prj[i] = prj[i].reshape((n_freqs, 16, -1))
    #t1 = time.time()
    tmp = predict(prj, x_tensors)
    #t2 = time.time()
    #print("time: ", t2-t1, tmp.shape)
    return tmp
