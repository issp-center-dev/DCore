from __future__ import print_function

from pytriqs.applications.dft.converters import *
from pytriqs.archive import *
from pytriqs.utility.h5diff import h5diff
import pytriqs.utility.mpi as mpi

from pytriqs.applications.pydmft.dmft_core import DMFTCoreSolver

#print("AAAAA")
params = {'aaa' : 0.1}

solver = DMFTCoreSolver('dmft', params)
#DMFTSolver
#Converter = Wannier90Converter(seedname='LaVO3-Pnma',hdf_filename='w90_convert.out.h5')
#
#Converter.convert_dft_input()
#
#if mpi.is_master_node():
    #h5diff("w90_convert.out.h5","w90_convert.ref.h5") 
