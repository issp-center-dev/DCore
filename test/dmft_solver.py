from __future__ import print_function

from pytriqs.applications.dft.converters import *
from pytriqs.archive import *
from pytriqs.utility.h5diff import h5diff
import pytriqs.utility.mpi as mpi

from pytriqs.applications.pydmft.standard import standard
from pytriqs.applications.pydmft.dmft_core import DMFTCoreSolver

#try:
    #import configparser
#except ImportError:
    #import ConfigParser as configparser

# Generate a HDF5 file
standard('stan.in')

#parser = configparser.ConfigParser()
#Make case sensitive
#parser.optionxform = str
#parser.read("dmft_solver.ini")

#params = {}
#params['system'] = {}

#print(parser._sections)
#
#params = {}
#for sect in parser.sections():
    #params[sect] = {}
    #for opt in parser.options(sect):
        #params[sect][opt] = parser.get(sect, opt)

params = {}
params['system'] = {}
params['system']['beta'] = 10.0
params['system']['n_iw'] = 2048
params['system']['n_tau'] = 10000
params['system']['dc_type'] = -1
params['system']['fix_mu'] = False

params['impurity_solver'] = {}
params['impurity_solver']['N_l'] = 50
params['impurity_solver']['name'] = 'TRIQS/hubbard-I'

params['control'] = {}
params['control']['sigma_mix'] = 0.5
params['control']['delta_mix'] = 0.5

solver = DMFTCoreSolver('test', params)

# FIXME: I want to use the Hubbard-I solver for tests!
solver.solve(max_step=1, output_file='test.out.h5', output_group='dmft_out', dry_run=True)

#DMFTSolver
#Converter = Wannier90Converter(seedname='LaVO3-Pnma',hdf_filename='w90_convert.out.h5')
#
#Converter.convert_dft_input()
#
#if mpi.is_master_node():
    #h5diff("w90_convert.out.h5","w90_convert.ref.h5") 
