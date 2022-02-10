from datetime import datetime
import numpy
from numpy.testing import assert_allclose
from itertools import product
from dcore._dispatcher import GfImFreq, BlockGf, HDFArchive

from .tools import float_to_complex_array, get_block_size
from .dcore_w90tool import Wannier90

def mk_hr_square_2x2(nspin, t, seedname):
    # spatial dimension: 1, 2, or 3
    dim = 2

    # Size of unit cell
    L = 2
    nrpts = 3**dim

    num_sites = L**dim

    Lx = 2 if dim>=1 else 1
    Ly = 2 if dim>=2 else 1
    Lz = 2 if dim>=3 else 1

    HamR_full = numpy.zeros([nrpts, nspin, num_sites, nspin, num_sites])
    irvec = numpy.empty((nrpts, 3), dtype=int)
    pos_in_unit_cell = numpy.array([ [i,j,k] for i, j, k in product(range(Lx), range(Ly), range(Lz))])

    ndgen = numpy.ones((nrpts,), dtype=int)

    Rx = 1 if dim>=1 else 0
    Ry = 1 if dim>=2 else 0
    Rz = 1 if dim>=3 else 0

    primitive_vecs = numpy.array([[Rx*2, 0, 0], [0, Ry*2, 0], [0, 0, Rz*2]])
    
    for ir, (X, Y, Z) in enumerate(product(range(-Rx,Rx+1), range(-Ry,Ry+1), range(-Rz,Rz+1))):
        irvec[ir, :] = (X, Y, Z)
    
        # <0 |H| R>
        for i_lsite, i_rsite in product(range(num_sites), repeat=2):
            pos_lsite = pos_in_unit_cell[i_lsite]
            pos_rsite = pos_in_unit_cell[i_rsite] + X * primitive_vecs[0] + Y * primitive_vecs[1] + Z * primitive_vecs[2]

            if numpy.sum((pos_lsite - pos_rsite)**2) == 1:
                for isp in range(nspin):
                    HamR_full[ir, isp, i_lsite, isp, i_rsite] = - t
    HamR_full = HamR_full.reshape((nrpts, nspin*num_sites, nspin*num_sites))

    print('Use the following line for input of DCore')
    print('corr_to_inequiv = ', ", ".join([str(numpy.sum(p)%2) for p in pos_in_unit_cell]))

    write_hr(F"""{seedname}_hr.dat""", irvec, HamR_full, ndgen)


def mk_hr_square(nf, t, seedname):
    """
    nf:
       Number of spin orbitals
    """
    # Size of unit cell
    dim = 2
    nrpts = 3**dim

    HamR = numpy.zeros([nrpts, nf, nf], dtype=numpy.complex128)
    irvec = numpy.empty((nrpts, 3), dtype=int)
    
    ndgen = numpy.ones((nrpts,), dtype=int)
    for ir, (X, Y) in enumerate(product(range(-1,2), range(-1,2))):
        irvec[ir, :] = (X, Y, 0)
        if numpy.sum(irvec[ir, :]**2) == 1:
            # <0 |H| R>
            HamR[ir, :, :] = - t * numpy.identity(nf)
    write_hr(F"""{seedname}_hr.dat""", irvec, HamR, ndgen)


# Write data to *_hr.dat
def write_hr(filename, irvec, HamR, ndgen):
    nrpts, norb = HamR.shape[0:2]
    with open(filename, 'w') as f:
        print(datetime.now(), file=f)
        print(norb, file=f)
        print(nrpts, file=f)
        for k in range(nrpts):
            print(ndgen[k], file=f, end=' ')
            if k % 15 == 14 or k == nrpts - 1:
                print('', file=f)
        for k in range(nrpts):
            for j, i in product(range(norb), repeat=2):
                print("{} {} {}  {} {}  {} {}".format(
                    irvec[k, 0], irvec[k, 1], irvec[k, 2],
                    i + 1, j + 1, HamR[k, i, j].real, HamR[k, i, j].imag), file=f)


def create_random_self_energy(block_names, dim_sh, beta, n_iw, noise):
    # Random self-energy
    Sigma_iw_sh = []
    for ish in range(len(dim_sh)):
        block_gfs = []
        for _ in block_names:
            s = GfImFreq(indices = numpy.arange(dim_sh[ish]).tolist(),
                beta = beta, n_points = n_iw, name = "sh{ish}")
            data_ = noise * numpy.random.randn(n_iw, dim_sh[ish], dim_sh[ish])
            data_ = data_ + data_.transpose((0, 2, 1))
            data_pn_ = numpy.empty((2*n_iw, dim_sh[ish], dim_sh[ish]))
            data_pn_[n_iw:, :, :] = data_
            data_pn_[0:n_iw, :, :] = data_[::-1, :, :]
            s.data[...] = data_pn_
            block_gfs.append(s)
        Sigma_iw_sh.append(BlockGf(name_list = block_names, block_list = block_gfs, make_copies = True))
    return Sigma_iw_sh


def block_to_numpy(blockgf_sh, bname, corr_to_inequiv):
    """
    Return a numpy array that contains the Green's function data for
    the target block.
    The data is ordered as crsh0, crsh1, .....

    blockgf_sh: list of BlockGf objects
       Green's function data for inequivalent shells
    bname: str
       Name of the target block
    corr_to_inequivalent: int 1D array
       Mappling from correlated shells to inequivalent shells
    """
    ncorr_sh = len(corr_to_inequiv)
    block_arr = [blockgf_sh[corr_to_inequiv[icorr]][bname].data for icorr in range(ncorr_sh)]
    nfreqs = block_arr[0].shape[0]
    totdim = numpy.sum([x.shape[1] for x in block_arr])
    bl_arr = numpy.zeros((nfreqs, totdim, totdim), dtype=numpy.complex128)
    offset = 0
    for icorr in range(ncorr_sh):
        offset_next = offset + block_arr[icorr].shape[1]
        bl_arr[:, offset:offset_next, offset:offset_next] = \
            block_arr[icorr]
        offset = offset_next
    return bl_arr

def blockgf_sh_to_numpy(blockgf_sh, corr_to_inequiv):
    """
    Return a numpy array that contains the Green's function data of
    all corrleted shells from given BlockGf objects defined on inequivalent shells.

    The spin orbitals are ordered as (crsh0, up), (crsh1, up), ..., (crsh0, down), (crsh1, down), ...

    blockgf_sh: list of BlockGf objects
       Green's function data for inequivalent shells
    corr_to_inequivalent: int 1D array
       Mappling from correlated shells to inequivalent shells
    """
    if blockgf_sh[0].n_blocks == 2:
        # spin_orbit = False
        bl_arr_u = block_to_numpy(blockgf_sh, 'up', corr_to_inequiv)
        bl_arr_d = block_to_numpy(blockgf_sh, 'down', corr_to_inequiv)
        nf = 2*bl_arr_u.shape[1]
        bl_arr_ud = numpy.zeros((bl_arr_u.shape[0], nf, nf), dtype=numpy.complex128)
        bl_arr_ud[:, 0:nf//2, 0:nf//2] = bl_arr_u
        bl_arr_ud[:, nf//2:, nf//2:] = bl_arr_d
    else:
        ncrsh = corr_to_inequiv.size
        nfreqs = blockgf_sh[0]['ud'].data.shape[0]
        num_so_inequiv_sh = numpy.array([get_block_size(bg['ud']) for bg in blockgf_sh])
        num_so_crsh = numpy.array([num_so_inequiv_sh[corr_to_inequiv[icrsh]] for icrsh in range(ncrsh)])
        num_orb = numpy.sum(num_so_crsh)//2
        bl_arr_ud = numpy.zeros((nfreqs, 2, num_orb, 2, num_orb), dtype=numpy.complex128)
        offset = 0
        for icrsh in range(ncrsh):
            offset_next = offset + num_so_crsh[icrsh]//2
            bl_arr_ud[:, :, offset:offset_next, :, offset:offset_next] = \
                blockgf_sh[corr_to_inequiv[icrsh]]['ud'].data.reshape((nfreqs, 2, num_so_crsh[icrsh]//2, 2, num_so_crsh[icrsh]//2))
            offset = offset_next
        bl_arr_ud = bl_arr_ud.reshape((bl_arr_ud.shape[0], 2*num_orb, 2*num_orb))
    return bl_arr_ud


def gk_from_w90(seedname, beta, nk1, nk2, nk3, Sigma_iw, smpl_freqs, corr_to_inequiv, mu=0.0):
    sigma_iw_numpy_dense = blockgf_sh_to_numpy(Sigma_iw, corr_to_inequiv)
    n_iw = sigma_iw_numpy_dense.shape[0]//2
    nf = sigma_iw_numpy_dense.shape[1]

    w90 = Wannier90(seedname+'_hr.dat')
    nwann = w90.Nwann

    # Compute one-particle Green's function on fermionic sampling frequencies
    nk = nk1 * nk2 * nk3
    givk = numpy.zeros((smpl_freqs.size, nf, nf, nk), dtype=numpy.complex128)
    ik1, ik2, ik3 = numpy.meshgrid(numpy.arange(nk1), numpy.arange(nk2), numpy.arange(nk3), indexing='ij')
    ik1, ik2, ik3 = ik1.ravel(), ik2.ravel(), ik3.ravel()
    iv = 1J * (2*smpl_freqs+1) * numpy.pi/beta
    hk = numpy.zeros((nk, nf, nf), dtype=numpy.complex128)
    if nwann == nf:
        for ik  in range(nk):
            hk[ik, :, :] = w90.get_Hk((ik1[ik]/nk1, ik2[ik]/nk2, ik3[ik]/nk3))
    elif 2*nwann == nf:
        for ik  in range(nk):
            hk[ik, 0:nwann, 0:nwann] = hk[ik, nwann:2*nwann, nwann:2*nwann] = \
                w90.get_Hk((ik1[ik]/nk1, ik2[ik]/nk2, ik3[ik]/nk3))
    for idx_w, iv_ in enumerate(iv):
        sigma = numpy.zeros((nf, nf), dtype=numpy.complex128)
        if abs(smpl_freqs[idx_w]) < n_iw:
            sigma[...] = sigma_iw_numpy_dense[smpl_freqs[idx_w] + n_iw,:,:]
            sigma -= mu * numpy.identity(nf)
        for ik in range(nk):
            givk[idx_w, :, :, ik] = numpy.linalg.inv(
                iv_ * numpy.identity(nf) - sigma - hk[ik, :, :])
    return givk


def gk_square(t, nf, beta, nk1, nk2, Sigma_iw, smpl_freqs, corr_to_inequiv, mu=0.0):
    sigma_iw_numpy_dense = blockgf_sh_to_numpy(Sigma_iw, corr_to_inequiv)
    n_iw = sigma_iw_numpy_dense.shape[0]//2

    # Compute one-particle Green's function on fermionic sampling frequencies
    nk = nk1 * nk2
    givk = numpy.zeros((smpl_freqs.size, nf, nf, nk), dtype=numpy.complex128)
    ik1, ik2 = numpy.meshgrid(numpy.arange(nk1), numpy.arange(nk2), indexing='ij')
    ik1, ik2 = ik1.ravel(), ik2.ravel()
    ek = -2*t*(numpy.cos(2*numpy.pi*ik1/nk1) + numpy.cos(2*numpy.pi*ik2/nk2))
    iv = 1J * (2*smpl_freqs+1) * numpy.pi/beta
    for idx_w, iv_ in enumerate(iv):
        sigma = numpy.zeros((nf, nf), dtype=numpy.complex128)
        if abs(smpl_freqs[idx_w]) < n_iw:
            sigma[...] = sigma_iw_numpy_dense[smpl_freqs[idx_w] + n_iw,:,:]
            sigma -= mu * numpy.identity(nf)
        for ik in range(nk1*nk2):
            givk[idx_w, :, :, ik] = numpy.linalg.inv(
                (iv_ - ek[ik]) * numpy.identity(nf) - sigma)
    return givk


def gk_tail(nf, beta, nk, smpl_freqs):
    """1/iw"""
    iw = 1J * (2 * smpl_freqs + 1) * numpy.pi/beta
    return numpy.einsum('k,w,ij->wijk', numpy.ones(nk), 1/iw, numpy.identity(nf))


def check_self_energy(seedname, Sigma_iw_sh, block_names, ite=1):
    n_inequiv_sh = len(Sigma_iw_sh)
    with HDFArchive(f"""{seedname}.out.h5""", 'r') as h:
        for ish in range(n_inequiv_sh):
            for bname in block_names:
                x = float_to_complex_array(h['dmft_out']['Sigma_iw'][f"""ite{ite}"""]
                    [f"""sh{ish}"""][bname]['data'])
                y = Sigma_iw_sh[ish][bname].data
                assert_allclose(x, y, rtol=0, atol=1e-10*numpy.abs(y).max())