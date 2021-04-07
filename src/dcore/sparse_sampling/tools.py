

import numpy
import copy

from irbasis_util.four_point_ph_view import FourPointPHView
from irbasis_util.four_point import FourPoint
from irbasis_util.tensor_regression import enable_MPI, fit, predict

from mpi4py import MPI 

enable_MPI() #for tensor_regression

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
is_master_node = comm.Get_rank() == 0

def construct_prj(phb, freqs_f_list):
    """

    Parameters
    ----------
    phb: FourPointPHView
        basis object
    freqs_list: list of (int, int)
        list of two fermionic frequencies
    x_tensors: list of ndarray
        Decomposed tensors

    Returns
    -------
    list of ndarray

    """
    assert isinstance(phb, FourPointPHView)
    assert isinstance(freqs_f_list, numpy.ndarray)
    assert freqs_f_list.dtype == int
    assert freqs_f_list.ndim == 2
    assert freqs_f_list.shape[1] == 2

    prj = phb.projector_to_matsubara_vec(freqs_f_list, decomposed_form=True)
    for i in range(len(prj)):
        prj[i] = prj[i].reshape((-1, 12, phb.Nl))
    return prj


def construct_prj_three_freqs(basis, freqs_list):
    """

    Parameters
    ----------
    phb: FourPoint
        basis object
    freqs_list: list of (int, int, int, int)
        list of four fermionic frequencies
    x_tensors: list of ndarray
        Decomposed tensors

    Returns
    -------
    list of ndarray

    """
    assert isinstance(basis, FourPoint)
    assert isinstance(freqs_list, numpy.ndarray)
    assert freqs_list.dtype == int
    assert freqs_list.ndim == 2
    assert freqs_list.shape[1] == 4

    prj = basis.projector_to_matsubara_vec(freqs_list, decomposed_form=True)
    for i in range(len(prj)):
        prj[i] = prj[i].reshape((-1, 16, basis.Nl))
    return prj

def perform_fit(projectors, xloc_local, D, niter, seed, alpha, rtol):
    num_o, num_freqs = xloc_local.shape
    #num_rep = projectors[0].shape[1] # Number of representations (maybe 12 or 16)
    #linear_dim = projectors[0].shape[-1] # dim of IR basis

    squared_norm_orb = comm.allreduce(numpy.sqrt(numpy.sum(numpy.abs(xloc_local)**2, axis=-1)).ravel())
    orb_idx = squared_norm_orb / numpy.amax(squared_norm_orb) > 1e-8
    if rank == 0:
        print("Num of active orbital components = ", numpy.sum(orb_idx))

    # (orbital, freq) => (freq, orbital)
    y = xloc_local[orb_idx, :].transpose((1,0))
    xs = fit(y, projectors, D, niter, rtol, alpha=alpha, verbose=1, random_init=True, comm=comm, seed=seed, nesterov=True)

    x_orb_full = numpy.zeros((D, num_o), dtype=complex)
    x_orb_full[:, orb_idx] = xs[-1]
    xs[-1] = x_orb_full

    return xs

def predict_xloc(prj, x_tensors):
    return predict(prj, x_tensors).transpose()


def interpolate(basis, freqs_list, x_tensors):
    """

    Parameters
    ----------
    basis: FourPointPHView or FourPoint
        basis object
    freqs_list: list of (int, int) or (int, int, int)
        list of two fermionic frequencies (FourPointPHView) or (******) (FourPoint)
    x_tensors: list of ndarray
        Decomposed tensors

    Returns
        2D array of (frequencies, spin-orbitals)
    -------


    """
    assert isinstance(basis, FourPointPHView) or isinstance(basis, FourPoint)
    assert isinstance(freqs_list, list)
    assert isinstance(x_tensors, list)

    prj = basis.projector_to_matsubara_vec(freqs_list, decomposed_form=True)
    return predict(prj, x_tensors)
