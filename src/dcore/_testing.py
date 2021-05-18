
from datetime import datetime
import numpy
from itertools import product

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