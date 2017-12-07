from __future__ import print_function

from pytriqs.applications.dft.converters import *
from pytriqs.archive import *
from pytriqs.utility.h5diff import h5diff
import pytriqs.utility.mpi as mpi

from pytriqs.applications.dcore.dcore_pre import dcore_pre
from pytriqs.applications.dcore.dmft_core import DMFTCoreSolver, create_parser

import tempfile

#try:
    #import configparser
#except ImportError:
    #import ConfigParser as configparser

seedname = "test1"

# Generate a HDF5 file
f = open('stan.in', 'w')
print("[model]", file=f)
print("t = 1.0", file=f)
print("U = 4.0", file=f)
print("seedname = "+seedname, file=f)
f.close()

dcore_pre('stan.in')

# Create default parameters
parser = create_parser()
tmp = tempfile.NamedTemporaryFile(mode='w+r')
parser.read(tmp.name)

params = parser.as_dict()
params['system']['beta'] = 10.0
params['system']['n_iw'] = 2048
params['system']['n_tau'] = 10000
params['system']['dc_type'] = -1
params['system']['fix_mu'] = False

params['impurity_solver']['n_l{int}'] = 50
params['impurity_solver']['name'] = 'TRIQS/hubbard-I'

params['control']['sigma_mix'] = 0.5
params['control']['delta_mix'] = 0.5

solver = DMFTCoreSolver(seedname, params)

# FIXME: I want to use the Hubbard-I solver for tests!
solver.solve(max_step=1, output_file=seedname+'.out.h5', output_group='dmft_out', dry_run=True)

