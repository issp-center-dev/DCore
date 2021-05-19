from datetime import datetime
import numpy
from itertools import product
from triqs.gf import GfImFreq, BlockGf

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