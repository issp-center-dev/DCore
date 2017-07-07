#!/usr/bin/env python
import sys
import numpy
from pytriqs.applications.dft.converters.hk_converter import *
from pytriqs.operators.util.U_matrix import *
#
# If input file is not specified ... 
#
if len(sys.argv) != 2:
    print "Usage:"
    print "$ ipytriqs standard.py input"
    sys.exit()
#
# Set Default value
#
t= 1.0
tp = 0.0
U = 0.0
J = 0.0
norb = 1
nk = [8, 1, 1]
lattice = "chain"
nelec = 1.0
seedname = "pydmft"
#
# Perse keywords and store
#
for line in open(sys.argv[1], 'r'):
    itemList = line.split('=')
    if itemList[0].strip() == 't':
        t = float(itemList[1])
    elif itemList[0].strip() == 't\'':
        tp = float(itemList[1])
    elif itemList[0].strip() == 'U':
        U = float(itemList[1])
    elif itemList[0].strip() == 'J':
        J = float(itemList[1])
    elif itemList[0].strip() == 'norb':
        norb = int(itemList[1])
    elif itemList[0].strip() == 'nk':
        nk[0] = int(itemList[1])
    elif itemList[0].strip() == 'lattice':
        lattice = itemList[1].strip()
    elif itemList[0].strip() == 'nelec':
        nelec = float(itemList[1])
    elif itemList[0].strip() == 'seedname':
        seedname = itemList[1].strip()
    else:
        print "Error ! Invalid keyword : ", itemList[0]
        sys.exit()
#
# Summary of input parameters
#
print "Parameter summary"
print "                 t = ", t
print "                t' = ", tp
print "                 U = ", U
print "                 J = ", J
print "Number of orbitals = ", norb
print "                nk = ", nk[0]
print "           Lattice = ", lattice.strip()
print "          seedname = ", seedname.strip()
print "             nelec = ", nelec
if lattice.strip() == 'chain':
    ndim = 1
elif lattice.strip() == 'square':
    ndim = 2
    nk[1] = nk[0]
elif lattice.strip() == 'cubic':
    ndim = 3
    nk[1] = nk[0]
    nk[2] = nk[0]
elif lattice.strip() == 'bethe':
    ndim = 0
else:
    print "Error ! Invalid lattice : ", lattice
    sys.exit()
nkBZ = nk[0] * nk[1] * nk[2]
print "         Dimension = ", ndim
print " Total number of k = ", nkBZ
#
# Write General-Hk formated file
#
f= open(seedname, 'w')
f.write(str(nkBZ)+"\n")
f.write(str(nelec)+"\n")
f.write("1\n")
f.write("1 1 0 " + str(norb)+"\n")
f.write("1\n")
f.write("1 1 0 " + str(norb) + " 0 0\n")
f.write("1 " + str(norb) +"\n")

kvec = [0.0, 0.0, 0.0]
for i0 in range(nk[0]):
    kvec[0] = 2.0 * numpy.pi * float(i0) / float(nk[0])
    for i1 in range(nk[1]):
      kvec[1] = 2.0 * numpy.pi * float(i1) / float(nk[1])
      for i2 in range(nk[2]):
           kvec[2] = 2.0 * numpy.pi * float(i2) / float(nk[2])
           ek = 2*t*(numpy.cos(kvec[0]) + numpy.cos(kvec[1]) + numpy.cos(kvec[2]))
           f.write(str(ek) + "\n")  #Real part
           f.write("0.0" + "\n") #Imaginary part    
f.close()
#
# Convert General-Hk to SumDFT-HDF5 format
#
Converter = HkConverter(filename = seedname)
Converter.convert_dft_input()
#
# Add U-matrix block
#
Umat, Upmat = U_matrix_kanamori(n_orb=norb, U_int=U, J_hund=J)
f = HDFArchive(seedname+'.h5','a')
if not ("pyDMFT" in f):
    f.create_group("pyDMFT")
f["pyDMFT"]["U_matrix"] = Umat
f["pyDMFT"]["Up_matrix"] = Upmat
#
# Finish
#
print "Done"
print ""