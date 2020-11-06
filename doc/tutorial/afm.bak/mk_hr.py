from __future__ import print_function
from datetime import datetime
from itertools import product
import numpy

# Size of unit cell
L = 2
nrpts = 3**3
t = 1
filename = 'cubic_hr.dat'

num_sites = L**3

HamR_full = numpy.zeros([nrpts, num_sites, num_sites])
irvec = numpy.empty((nrpts, 3), dtype=int)
pos_in_unit_cell = numpy.array([ [i,j,k] for i, j, k in product(range(2), repeat=3)])

ndgen = numpy.ones((nrpts,), dtype=int)

primitive_vecs = numpy.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])

for ir, (X, Y, Z) in enumerate(product(range(-1,2), repeat=3)):
    irvec[ir, :] = (X, Y, Z)

    # <0 |H| R>
    for i_lsite, i_rsite in product(range(num_sites), repeat=2):
        pos_lsite = pos_in_unit_cell[i_lsite]
        pos_rsite = pos_in_unit_cell[i_rsite] + X * primitive_vecs[0] + Y * primitive_vecs[1] + Z * primitive_vecs[2]

        if numpy.sum((pos_lsite - pos_rsite)**2) == 1:
            HamR_full[ir, i_lsite, i_rsite] = - t

print('Use the following line for input of DCore')
print('corr_to_inequiv = ', ", ".join([str(numpy.sum(p)%2) for p in pos_in_unit_cell]))

# Write data to *_hr.dat
with file(filename, 'w') as f:
    print(datetime.now(), file=f)
    print(num_sites, file=f)
    print(nrpts, file=f)
    for k in range(nrpts):
        print(ndgen[k], file=f, end=' ')
        if k % 15 == 14 or k == nrpts - 1:
            print('', file=f)
    
    for k in range(nrpts):
        for j, i in product(range(num_sites), repeat=2):
            print("{} {} {}  {} {}  {} {}".format(
                irvec[k, 0], irvec[k, 1], irvec[k, 2],
                i + 1, j + 1, HamR_full[k, i, j], 0.0), file=f)
