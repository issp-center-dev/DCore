#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy
from pytriqs.operators.util.U_matrix import U_matrix_kanamori
from pytriqs.archive.hdf_archive import HDFArchive
from pytriqs.applications.dft.converters.wannier90_converter import Wannier90Converter
from pytriqs.applications.dft.converters.hk_converter import HkConverter

def standard(filename):

    print("Reading {0} ...\n".format(filename))
    #
    # Set Default value
    #
    t= 1.0
    tp = 0.0
    U = 0.0
    J = 0.0
    model = "single"
    nk = 8
    nk0 = 0
    nk1 = 0
    nk2 = 0
    ncor = 1
    lattice = "chain"
    nelec = 1.0
    seedname = "pydmft"
    l = [0 for i in range(1000)]
    norb = [1 for i in range(1000)]
    equiv = [-1 for i in range(1000)]
    #
    # Perse keywords and store
    #
    for line in open(filename, 'r'):
        itemList = line.split('=')
        keyList = itemList[0].split('-')
        if itemList[0].strip() == 't':
            t = float(itemList[1])
        elif itemList[0].strip() == 't\'':
            tp = float(itemList[1])
        elif itemList[0].strip() == 'U':
            U = float(itemList[1])
        elif itemList[0].strip() == 'J':
            J = float(itemList[1])
        elif itemList[0].strip() == 'nk':
            nk = int(itemList[1])
        elif itemList[0].strip() == 'nk0':
            nk0 = int(itemList[1])
        elif itemList[0].strip() == 'nk1':
            nk1 = int(itemList[1])
        elif itemList[0].strip() == 'nk2':
            nk2 = int(itemList[1])
        elif itemList[0].strip() == 'ncor':
            ncor = int(itemList[1])
        elif itemList[0].strip() == 'lattice':
            lattice = itemList[1].strip()
        elif itemList[0].strip() == 'model':
            model = itemList[1].strip()
        elif itemList[0].strip() == 'nelec':
            nelec = float(itemList[1])
        elif itemList[0].strip() == 'seedname':
            seedname = itemList[1].strip()
        elif keyList[0].strip() == 'l':
            l[int(keyList[1])] = int(itemList[1])
        elif keyList[0].strip() == 'norb':
            norb[int(keyList[1])] = int(itemList[1])
        elif keyList[0].strip() == 'equiv':
            equiv[int(keyList[1])] = int(itemList[1])
        elif itemList[0].strip() == '':
            pass
        else:
            print("Error ! Invalid keyword : ", itemList[0])
            sys.exit()
    #
    # Summary of input parameters
    #
    print("Parameter summary")
    print("           Lattice = {0}".format(lattice.strip()))
    print("          seedname = {0}".format(seedname.strip()))
    print("                 U = {0}".format(U))
    print("                 J = {0}".format(J))
    print("             nelec = {0}".format(nelec))
    #
    # Lattice
    #
    if lattice.strip() == 'wannier90':
        if nk0 == 0: nk0 = nk
        if nk1 == 0: nk1 = nk
        if nk2 == 0: nk2 = nk
        print("                nk0 = {0}".format(nk0))
        print("                nk1 = {0}".format(nk1))
        print("                nk2 = {0}".format(nk2))
        print("               ncor = {0}".format(ncor))
        for i in range(ncor):
            if equiv[i] == -1: equiv[i] = i
            print("     l[{0}], norb[{0}], equiv[{0}] = {1}, {2}, {3}".format(i,l[i],norb[i],equiv[i]))
        #
        f = open(seedname+'.inp', 'w')
        print("0 {0} {1} {2}".format(nk0,nk1,nk2), file=f)
        print("{0}".format(nelec), file=f)
        print("{0}".format(ncor), file=f)
        for i in range(ncor):
            print("{0} {1} {2} {3} 0 0".format(i,equiv[i],l[i],norb[i]), file=f)
        f.close()
    else:
        print("                 t = {0}".format(t))
        print("                t' = {0}".format(tp))
        print("                nk = {0}".format(nk))
        print("             Model = {0}".format(model.strip()))
        weights_in_file = False
        if lattice.strip() == 'chain':
            nkBZ = nk
        elif lattice.strip() == 'square':
            nkBZ = nk**2
        elif lattice.strip() == 'cubic':
            nkBZ = nk**3
        elif lattice.strip() == 'bethe':
            nkBZ = nk
            weights_in_file = True
        else:
            print("Error ! Invalid lattice : ", lattice)
            sys.exit()
        print(" Total number of k =", str(nkBZ))
        #
        # Model
        #
        if model.strip() == 'single':
            l[0] = 0
            norb[0] = 1
        elif model.strip() == 'eg':
            #FIXME: l=2 does not make sense. l=2 assumes norb=5 (full d-shell) in generating Coulomb tensor.
            #What is the proper way to generate Colomb tensor for eg?
            l[0] = 2
            norb[0] = 2
        elif model.strip() == 't2g':
            l[0] = 1
            norb[0] = 3
        elif model.strip() == 'full-d':
            l[0] = 2
            norb[0] = 5
        else:
            print("Error ! Invalid lattice : {0}".format(model.strip()))
            sys.exit()
        #
        # Write General-Hk formated file
        #
        f = open(seedname, 'w')
        print("{0}".format(nkBZ), file=f)
        print("{0}".format(nelec), file=f)
        print("1", file=f)
        print("0 0 {0} {1}".format(l[0], norb[0]), file=f)
        print("1", file=f)
        print("0 0 {0} {1} 0 0".format(l[0], norb[0]), file=f)
        print("1 {0}".format(norb[0]), file=f)
        #
        # Energy band
        #
        if lattice.strip() == 'bethe':
            #
            # If Bethe lattice, set k-weight manually to generate semi-circular DOS
            #
            for i0 in range(nk):
                ek = float(2*i0 + 1 - nk) / float(nk)
                wk = numpy.sqrt(1.0 - ek**2)
                print("{0}".format(wk), file=f)
            for i0 in range(nk):
                ek = t * float(2*i0 + 1 - nk) / float(nk)
                for iorb in range(norb[0]):
                    for jorb in range(norb[0]):
                        if iorb == jorb:
                            print("{0}".format(ek), file=f) #Real part
                        else:
                            print("0.0", file=f) #Real part
                for iorb in range(norb[0]*norb[0]): print("0.0", file=f) #Imaginary part
        elif lattice.strip() == 'chain':
            kvec = [0.0, 0.0, 0.0]
            for i0 in range(nk):
                kvec[0] = 2.0 * numpy.pi * float(i0) / float(nk)
                ek = 2.0*t*numpy.cos(kvec[0]) + 2*tp*numpy.cos(2.0*kvec[0])
                for iorb in range(norb[0]):
                    for jorb in range(norb[0]):
                        if iorb == jorb:
                            print("{0}".format(ek), file=f) #Real part
                        else:
                            print("0.0", file=f) #Real part
                for iorb in range(norb[0]*norb[0]): print("0.0", file=f) #Imaginary part
        elif lattice.strip() == 'square':
            kvec = [0.0, 0.0, 0.0]
            for i0 in range(nk):
                kvec[0] = 2.0 * numpy.pi * float(i0) / float(nk)
                for i1 in range(nk):
                    kvec[1] = 2.0 * numpy.pi * float(i1) / float(nk)
                    ek = 2.0*t*(numpy.cos(kvec[0]) +  numpy.cos(kvec[1])) \
                       + 2.0*tp*(numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1]))
                    for iorb in range(norb[0]):
                        for jorb in range(norb[0]):
                            if iorb == jorb:
                                print("{0}".format(ek), file=f) #Real part
                            else:
                                print("0.0", file=f) #Real part
                    for iorb in range(norb[0]*norb[0]): print("0.0", file=f) #Imaginary part
        elif lattice.strip() == 'cubic':
            kvec = [0.0, 0.0, 0.0]
            for i0 in range(nk):
                kvec[0] = 2.0 * numpy.pi * float(i0) / float(nk)
                for i1 in range(nk):
                    kvec[1] = 2.0 * numpy.pi * float(i1) / float(nk)
                    for i2 in range(nk):
                        kvec[2] = 2.0 * numpy.pi * float(i2) / float(nk)
                        ek = 2*t*(numpy.cos(kvec[0]) +  numpy.cos(kvec[1]) + numpy.cos(kvec[2])) \
                           + 2*tp*( numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1]) \
                                  + numpy.cos(kvec[1] + kvec[2]) + numpy.cos(kvec[1] - kvec[2]) \
                                  + numpy.cos(kvec[2] + kvec[0]) + numpy.cos(kvec[2] - kvec[0]) )
                        for iorb in range(norb[0]):
                            for jorb in range(norb[0]):
                                if iorb == jorb:
                                    print("{0}".format(ek), file=f) #Real part
                                else:
                                    print("0.0", file=f) #Real part
                        for iorb in range(norb[0]*norb[0]): print("0.0", file=f) #Imaginary part
        f.close()
    #
    # Convert General-Hk to SumDFT-HDF5 format
    #
    if lattice.strip() == 'wannier90':
        Converter = Wannier90Converter(seedname = seedname)
        Converter.convert_dft_input()
    else:
        Converter = HkConverter(filename = seedname)
        Converter.convert_dft_input(weights_in_file=weights_in_file)
    #
    # Add U-matrix block (Tentative)
    # ####  The format of this block is not fixed  #### 
    #
    #Umat, Upmat = U_matrix_kanamori(n_orb=norb[0], U_int=U, J_hund=J)
    f = HDFArchive(seedname+'.h5','a')
    if not ("pyDMFT" in f):
        f.create_group("pyDMFT")
    #f["pyDMFT"]["U_matrix"] = Umat
    #f["pyDMFT"]["Up_matrix"] = Upmat
    f["pyDMFT"]["U_int"] = U
    f["pyDMFT"]["J_hund"] = J
    #
    # Finish
    #
    print("Done\n")

if __name__ == '__main__':

    #
    # If input file is not specified ... 
    #
    if len(sys.argv) != 2:
        print("Usage:")
        print("$ ipytriqs standard.py input")
        sys.exit()
    standard(sys.argv[1])
