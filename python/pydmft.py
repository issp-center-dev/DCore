#!/usr/bin/env python
import sys
from pytriqs.applications.dft.sumk_dft import *
#
# If input file is not specified ... 
#
if len(sys.argv) != 2:
    print "Usage:"
    print "$ ipytriqs pydmft.py input"
    sys.exit()
#
# Set Default value
#
seedname = "pydmft"
max_step = 10
#
# Perse keywords and store
#
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
from pytriqs.applications.impurity_solvers.hubbard_I import *
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