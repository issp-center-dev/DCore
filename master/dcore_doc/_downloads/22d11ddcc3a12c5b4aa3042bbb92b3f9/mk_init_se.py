from __future__ import print_function

import numpy
from itertools import product

norb = 1
mag = 1.0

with open('init_se_up.txt', 'w') as f:
    for i, (isp, iorb) in enumerate(product(range(2), range(norb))):
        print(isp, iorb, iorb, mag*(-1)**isp, 0.0, file=f)

with open('init_se_down.txt', 'w') as f:
    for i, (isp, iorb) in enumerate(product(range(2), range(norb))):
        print(isp, iorb, iorb, -mag*(-1)**isp, 0.0, file=f)
