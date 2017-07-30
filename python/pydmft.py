#!/usr/bin/env python
from __future__ import print_function
import sys, os
from pytriqs.applications.dft.sumk_dft import *

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

#
# If input file is not specified ... 
#
if len(sys.argv) != 2:
    print("Usage:")
    print("$ ipytriqs pydmft.py input")
    sys.exit()
#
# Set Default value
#
seedname = "pydmft"
max_step = 10
sigma_mix = 0.8                  # Mixing factor of Sigma after solution of the AIM
delta_mix = 0.8                  # Mixing factor of Delta as input for the AIM
dc_type = -1                     # DC type: -1 None, 0 FLL, 1 Held, 2 AMF
use_blocks = False               # use bloc structure from DFT input. If use_blocks == True, all Green's functions will be spin diagonal.

h_field = 0.0
prec_mu = 0.0001

n_iw=2049                        # Number of positive Matsubara frequencies
n_tau=20001                      # Number of tau points
n_l=50                           # Number of Legendre polynomials. FIXME: This should be an option of an impurity solver

solver = "TRIQS/cthyb"
#solver = "ALPS/cthyb"

#
# Parse keywords and store
#
ini_file = sys.argv[1]
ini = configparser.SafeConfigParser()
ini.optionxform = str #Distinguish uppercase letters and lowercase letters in option names
if os.path.exists(ini_file):
    ini.read(ini_file)
else:
    print("Input file not found: %s"%ini_file)
    sys.exit(1)

for line in open(sys.argv[1], 'r'):
    itemList = line.split('=')
    if itemList[0].strip() == 'seedname':
        seedname = itemList[1].strip()
    elif itemList[0].strip() == 'max_step':
        max_step = int(itemList[1])
    else:
        print "Error ! Invalid keyword : ", itemList[0]
        sys.exit()
#
# Summary of input parameters
#
print "Parameter summary"
print "          seedname = ", seedname.strip()
print "        Max. steps = ", max_step
#
# Open file for 1-body and U-matrix
#
SK = SumkDFT(hdf_file=seedname+'.h5')
U_file = HDFArchive(seedname+'.h5','r')
Umat = U_file["pyDMFT"]["U_matrix"]
Upmat = U_file["pyDMFT"]["Up_matrix"]
#
# Construct Green's function and impurity Hamiltonian
#
gf_struct = SK.gf_struct_solver[0]
h_int = h_int_density(spin_names, orb_names, map_operator_structure=SK.sumk_to_solver[0], U=Umat, Uprime=Upmat, H_dump="H.txt")
#
# Setup solver
#
from hubbard_solver_l0 import *
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=n_iw, n_tau=n_tau, n_l=n_l)
#
# DMFT roop
#
for iteration_number in range(max_step):
    #
    # Construct Local Green's function
    #
    #
    # Solve the impurity problem:
    #
    S.solve(h_int=h_int, **p)
    #
    # Convergence check
    #
    #
    # Update Green's function and Self energy
    #
#
# Finish
#
print "Done"
print ""