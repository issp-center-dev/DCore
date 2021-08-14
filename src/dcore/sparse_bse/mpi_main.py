import argparse
import sys
import numpy

import irbasis_x
from irbasis_x import bse as ir_bse
from irbasis_x import bse_dmft as bse_dmft
from mpi4py import MPI

from h5 import HDFArchive

from triqs.utility import mpi

from dcore.mpi import get_slice
from dcore.tools import float_to_complex_array, complex_to_float_array
from dcore.irbasis_util import fit_iw, construct_basis

def compute_glk(gk_file):
    comm = MPI.COMM_WORLD
    gkw, smpl_freqs, Lambda, cutoff, beta = None, None, None, None, None
    with HDFArchive(gk_file,  'r') as h:
        gkw = float_to_complex_array(h['data'])
        smpl_freqs = h['smpl_freqs']
        Lambda = h['Lambda_IR']
        cutoff = h['cutoff_IR']
        Lambda = h['Lambda_IR']
        cutoff = h['cutoff_IR']
        beta = h['beta']
    gkw = comm.bcast(gkw)
    smpl_freqs = comm.bcast(smpl_freqs)
    Lambda = comm.bcast(Lambda)
    cutoff = comm.bcast(cutoff)

    basis_f = construct_basis('F', beta, Lambda, cutoff)
    smpl_freqs_ = basis_f.sampling_points_matsubara(basis_f.dim()-1)
    assert numpy.array_equal(smpl_freqs, smpl_freqs_)

    nspin, norb = gkw.shape[1], gkw.shape[2]
    assert nspin == 2
    nf = 2*norb
    nk1, nk2, nk3 = gkw.shape[-3:]
    # (freq, spin, orb, spin, orb, nk1, nk2, nk3)
    gkw = gkw.reshape((gkw.shape[0], 2, nf//2, 2, nf//2, nk1, nk2, nk3))
    # (freq, orb, spin, orb, spin, nk1, nk2, nk3)
    gkw = gkw.transpose((0, 2,1, 4,3, 5, 6, 7))
    # (nl, orb, spin, orb, spin, nk1, nk2, nk3)
    glk = fit_iw(basis_f, gkw, axis=0)
    return glk, basis_f

def compute_Floc_from_G2loc(g2loc_file, gloc_l, basis_f):
    with HDFArchive(g2loc_file, 'r') as h:
        nsh = h['n_inequiv_sh']
        nso_sh = h['nso_sh']
        Lambda_IR = h['Lambda_IR']
        cutoff_IR = h['cutoff_IR']
        corr_to_inequiv = h['corr_to_inequiv']
        g2loc_sh = []
        wsample_ph = h['wsample_ph']
        wb_sample = h['wb_sample']
        for ish in range(nsh): # inequivalent shell
            g2loc_sh.append(
                float_to_complex_array(h['data'][f'sh{ish}'])
            )

    # Read Floc for each inequivalent shell
    n_corr_sh = corr_to_inequiv.size
    nsh = len(g2loc_sh)
    nso_corr_sh = numpy.array([nso_sh[corr_to_inequiv[i]] for i in range(n_corr_sh)])
    nso = numpy.sum(nso_corr_sh)
    nfreqs = g2loc_sh[0].shape[0]
    norb = nso//2
    nl = gloc_l.shape[0]

    # Mappling from an inequivalent shell to one of the corresponding correlated shells
    inequiv_to_corr = numpy.array([
        numpy.where(corr_to_inequiv==ish)[0][0] for ish in range(nsh)
    ])
    assert all([corr_to_inequiv[inequiv_to_corr[ish]] == ish for ish in range(nsh)])

    # Compute full vertex for each inequivalent shell
    basis = irbasis_x.FourPointBasis(Lambda_IR, basis_f.beta, cutoff_IR, vertex=True)
    Floc_sh = []
    assert gloc_l.shape == (nl, norb, 2, norb, 2)
    offset_corr_sh = numpy.hstack((0, numpy.cumsum(nso_corr_sh//2)))
    for ish in range(nsh):
        icrsh = inequiv_to_corr[ish]
        nso_sh = nso_corr_sh[icrsh]
        norb_sh = nso_corr_sh[icrsh]//2
        s_ = slice(offset_corr_sh[icrsh], offset_corr_sh[icrsh+1])

        ## Two-particle gf projected on the shell
        g2loc_ = g2loc_sh[ish]
        assert g2loc_.shape == (nfreqs,) + (2, norb_sh) * 4
        g2loc_ = g2loc_.transpose((0, 2,1, 4,3, 6,5, 8,7)) # (freq,) + (orb, spin)*4
        g2loc_ = g2loc_.reshape((nfreqs, nso_sh, nso_sh, nso_sh, nso_sh))

        ## One-particle gf projected on the shell
        gloc_sh_l = gloc_l[:, s_, :, s_, :] # (freq, orb, spin, orb, spin)
        gloc_sh_l = gloc_sh_l.reshape((nl, nso_sh, nso_sh)) # (freq, orb * spin, orb * spin)

        Floc_sh.append(
            ir_bse.compute_F_ph(
                    g2loc_sh[ish].reshape((nfreqs, nso_sh, nso_sh, nso_sh, nso_sh)),
                    wsample_ph, gloc_sh_l, basis_f
                )
        )

    return Floc_sh, wb_sample, wsample_ph, basis, corr_to_inequiv


def run(input_file, gk_file, g2loc_file, output_file):
    comm = MPI.COMM_WORLD

    # Input file
    with HDFArchive(input_file,  'r') as h:
        qsample = h['qsample']

    # Distribute q points over MPI processes
    slice_local = get_slice(qsample[0].size)
    qsample_local = \
        qsample[0][slice_local], qsample[1][slice_local], qsample[2][slice_local]

    # Fit one-particle gf with IR
    # glk: (nl, orb, spin, orb, spin, nk1, nk2, nk3)
    glk, basis_f = compute_glk(gk_file)
    nl = glk.shape[0]
    norb = glk.shape[1]
    nk1, nk2, nk3 = glk.shape[-3:]
    nf = 2*norb

    # Compute Floc
    # gloc_l: (nl, orb, spin, orb, spin)
    gloc_l = numpy.mean(glk, axis=(-3,-2,-1))
    Floc_sh, wb_sample, wsample_Floc_ph, basis, corr_to_inequiv = \
        compute_Floc_from_G2loc(g2loc_file, gloc_l, basis_f)

    nq = qsample[0].size
    num_wb = wb_sample.size
    chi = numpy.zeros((num_wb, nq, nf, nf, nf, nf), dtype=numpy.complex128)

    # Reshape one-particle gf
    glk = glk.reshape((nl, nf, nf, nk1, nk2, nk3))
    gloc_l = gloc_l.reshape((nl, nf, nf))

    offset = 0
    for idx_wb, wb in enumerate(wb_sample):
        # Solver
        solver = bse_dmft.SparseSolver(basis, 2*wb)
        nsmpl_freqs_  = solver.wsample_Floc[0].size
        offset_next = offset + nsmpl_freqs_
        for i in range(3):
            assert numpy.array_equal(
                wsample_Floc_ph[i][offset:offset_next], solver.wsample_Floc[i])

        # Compute X0q and X0loc
        X0q_ph = bse_dmft.compute_X0q_ph(glk, basis_f, solver.wsample_X0, qsample_local)
        X0loc = bse_dmft.compute_X0loc_ph(glk, basis_f, solver.wsample_X0)

        # Project Floc_sh to wb
        Floc_sh_wb = [x[offset:offset_next,...] for x in Floc_sh]

        # chi_local: (nf, nf, nf, nf, nq_local)
        chi_local_ = solver.solve(X0q_ph, X0loc, Floc_sh_wb, corr_to_inequiv)
        chi[idx_wb, slice_local, ...] = numpy.moveaxis(chi_local_, -1, 0)

        offset = offset_next

    chi = comm.allreduce(chi)
    chi = chi.reshape((num_wb, nq) + (nf//2, 2) * 4)
    chi = chi.transpose((0,1, 3,2, 5,4, 7,6, 9,8)) # (num_wb,nq) + (2,norb)*4

    # Save results into a HDF5 file
    if mpi.is_master_node():
        with HDFArchive(output_file, 'w') as h:
            h['nk'] = numpy.array([nk1, nk2, nk3])
            h['qsample'] = qsample
            h['wb_sample'] = wb_sample
            h['chi'] = complex_to_float_array(chi)
            h['wsample_Floc_ph'] = wsample_Floc_ph
            h.create_group('Floc_sh')
            for ish, Floc_sh_ in enumerate(Floc_sh):
                h['Floc_sh'][f'sh{ish}'] = complex_to_float_array(Floc_sh_)
    return


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(
            description='MPI program for solving BSE')
        parser.add_argument('input_file')
        parser.add_argument('gk_file')
        parser.add_argument('g2loc_file')
        parser.add_argument('output_file')
        args = parser.parse_args()

        run(args.input_file, args.gk_file, args.g2loc_file, args.output_file)


    except Exception as e:
        print("Unexpected error:", e)
        import traceback
        traceback.print_exc()
        sys.exit(1)
