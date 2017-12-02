#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import numpy
import argparse
import re
from pytriqs.archive.hdf_archive import HDFArchive
from pytriqs.applications.dft.converters.wannier90_converter import Wannier90Converter
from pytriqs.applications.dft.converters.hk_converter import HkConverter
from dmft_core import create_parser
import pytriqs.utility.mpi as mpi

#from .typed_parser import TypedParser
from typed_parser import TypedParser


def __print_paramter(p, param_name):
    print(param_name + " = " + str(p[param_name]))


def __generate_wannier90_model(params, l, norb, equiv, f):
    """
    Compute hopping etc. for A(k,w) of Wannier90

    Parameters
    ----------
    params : dictionary
        Input parameters
    l : integer array
        Angular momentum at each correlation shell
    norb : integer array
        Number of orbitals at each correlation shell
    equiv : integer array
        Equivalence of correlation shell
    f : File stream
        file
    """
    nk = params["system"]["nk"]
    nk0 = params["system"]["nk0"]
    nk1 = params["system"]["nk1"]
    nk2 = params["system"]["nk2"]
    ncor = params["model"]["ncor"]

    if nk0 == 0: nk0 = nk
    if nk1 == 0: nk1 = nk
    if nk2 == 0: nk2 = nk
    print("\n    nk0 = {0}".format(nk0))
    print("    nk1 = {0}".format(nk1))
    print("    nk2 = {0}".format(nk2))
    print("    ncor = {0}".format(ncor))
    for i in range(ncor):
        assert equiv[i] >= 0
        print("    l[{0}], norb[{0}], equiv[{0}] = {1}, {2}, {3}".format(i,l[i],norb[i],equiv[i]))
    print("")
    #
    # Generate file for Wannier90-Converter
    #
    print("0 {0} {1} {2}".format(nk0,nk1,nk2), file=f)
    print(params["model"]["nelec"], file=f)
    print(ncor, file=f)
    for i in range(ncor):
        print("{0} {1} {2} {3} 0 0".format(i,equiv[i],l[i],norb[i]), file=f)


def __generate_lattice_model(params, l, norb, equiv, f):
    """
    Compute hopping etc. for A(k,w) of preset models

    Parameters
    ----------
    params : dictionary
            Input parameters
    Returns
    -------
    weights_in_file: Bool
        Weights for the bethe lattice
    """
    weights_in_file = False
    if params["model"]["lattice"] == 'chain':
        nkBZ = params["system"]["nk"]
    elif params["model"]["lattice"] == 'square':
        nkBZ = params["system"]["nk"]**2
    elif params["model"]["lattice"] == 'cubic':
        nkBZ = params["system"]["nk"]**3
    elif params["model"]["lattice"] == 'bethe':
        nkBZ = params["system"]["nk"]
        weights_in_file = True
    else:
        print("Error ! Invalid lattice : ", params["model"]["lattice"])
        sys.exit()
    print("\n    Total number of k =", str(nkBZ))
    #
    # Model
    #
    if params["model"]["orbital_model"] == 'single':
        l[0] = 0
        norb[0] = 1
    elif params["model"]["orbital_model"] == 'eg':
        #FIXME: l=2 does not make sense. l=2 assumes norb=5 (full d-shell) in generating Coulomb tensor.
        #What is the proper way to generate Coulomb tensor for eg?
        l[0] = 2
        norb[0] = 2
    elif params["model"]["orbital_model"] == 't2g':
        l[0] = 1
        norb[0] = 3
    elif params["model"]["orbital_model"] == 'full-d':
        l[0] = 2
        norb[0] = 5
    else:
        print("Error ! Invalid lattice : ", params["model"]["orbital_model"])
        sys.exit()
    #
    # Write General-Hk formatted file
    #
    if mpi.is_master_node():
        print(nkBZ, file=f)
        print(params["model"]["nelec"], file=f)
        print("1", file=f)
        print("0 0 {0} {1}".format(l[0], norb[0]), file=f)
        print("1", file=f)
        print("0 0 {0} {1} 0 0".format(l[0], norb[0]), file=f)
        print("1 {0}".format(norb[0]), file=f)

    t = params["model"]["t"]
    tp = params["model"]["t'"]
    nk = params["system"]["nk"]

    #
    # Energy band
    #
    if params["model"]["lattice"] == 'bethe':
        #
        # If Bethe lattice, set k-weight manually to generate semi-circular DOS
        #
        for i0 in range(nk):
            ek = float(2*i0 + 1 - nk) / float(nk)
            wk = numpy.sqrt(1.0 - ek**2)
            print("{0}".format(wk), file=f)
        for i0 in range(nk):
            ek = 2.0 * t * float(2*i0 + 1 - nk) / float(nk)
            for iorb in range(norb[0]):
                for jorb in range(norb[0]):
                    if iorb == jorb:
                        print("{0}".format(ek), file=f) #Real part
                    else:
                        print("0.0", file=f) #Real part
            for iorb in range(norb[0]*norb[0]): print("0.0", file=f) #Imaginary part
    elif params["model"]["lattice"] == 'chain':
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
    elif params["model"]["lattice"] == 'square':
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
    elif params["model"]["lattice"] == 'cubic':
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
    return weights_in_file


def pydmft_pre(filename):
    """
    Main routine for the pre-processing tool

    Parameters
    ----------
    filename : string
        Input-file name
    """
    if mpi.is_master_node(): print("\n  @ Reading {0} ...".format(filename))
    #
    # Construct a parser with default values
    #
    parser = create_parser()
    #
    # Parse keywords and store
    #
    parser.read(filename)
    p = parser.as_dict()
    #
    # cshell=(l, norb, equiv) or (l, norb)
    #
    cshell_list=re.findall(r'\(\s*\d+,\s*\d+,*\s*\d*\)', p["model"]["cshell"])
    l = [0]*p["model"]['ncor']
    norb = [1]*p["model"]['ncor']
    equiv = [-1]*p["model"]['ncor']
    try:
        equiv_str_list = []
        equiv_index = 0
        for  i, _list  in enumerate(cshell_list):
            _cshell = filter(lambda w: len(w) > 0, re.split(r'[\(\s*\,\s*,*\s*\)]', _list))
            l[i] = int(_cshell[0])
            norb[i] = int(_cshell[1])
            if len(_cshell) == 3:
                if _cshell[2] in equiv_str_list:
                    # Defined before
                    equiv[i] = equiv_str_list.index(_cshell[2])
                else:
                    # New one
                    equiv_str_list.add(_cshell[2])
                    equiv[i] = equiv_index
                    ++ equiv_index
            else:
                equiv[i] = equiv_index
                ++ equiv_index
    except:
        raise RuntimeError("Error ! Format of cshell is wrong.")
    #
    # Summary of input parameters
    #
    if mpi.is_master_node():
        print("\n  @ Parameter summary")
        print("\n    [model] block")
        for k, v in p["model"].items():
            print("      {0} = {1}".format(k, v))
        print("\n    [system] block")
        for k, v in p["system"].items():
            print("      {0} = {1}".format(k, v))
    #
    # Lattice
    #
    if mpi.is_master_node():
        print("\n  @ Generate model-HDF5 file")
    seedname = p["model"]["seedname"]
    if p["model"]["lattice"] == 'wannier90':
        if mpi.is_master_node():
            with open(seedname+'.inp', 'w') as f:
                __generate_wannier90_model(p, l, norb, equiv, f)
        # Convert General-Hk to SumDFT-HDF5 format
        Converter = Wannier90Converter(seedname = seedname)
        Converter.convert_dft_input()
    else:
        if mpi.is_master_node():
            with open(seedname+'.inp', 'w') as f:
                weights_in_file = __generate_lattice_model(p, l, norb, equiv, f)
        # Convert General-Hk to SumDFT-HDF5 format
        Converter = HkConverter(filename = seedname + ".inp", hdf_filename=seedname+".h5")
        Converter.convert_dft_input(weights_in_file=weights_in_file)
    #
    # Add U-matrix block (Tentative)
    # ####  The format of this block is not fixed  ####
    #
    if mpi.is_master_node():
        print("\n  @ Write the information of interactions")
        f = HDFArchive(seedname+'.h5','a')
        if not ("pyDMFT" in f):
            f.create_group("pyDMFT")
        f["pyDMFT"]["U_int"] = p["model"]["U"]
        f["pyDMFT"]["J_hund"] = p["model"]["J"]
        print("\n    Wrote to {0}".format(seedname+'.h5'))
    #
    # Finish
    #
    if mpi.is_master_node(): print("\n  Done\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(\
        prog='pydmft_pre.py',\
        description='pre script for pydmft.',\
        epilog='end',\
        usage = '$ pydmft_pre input',\
        add_help= True)
    parser.add_argument('path_input_file', \
                        action = 'store',\
                        default= None,    \
                        type=str, \
                        help = "input file name."
    )
    
    args=parser.parse_args()
    if(os.path.isfile(args.path_input_file) is False):
        print("Input file is not exist.")
        sys.exit()
    pydmft_pre(args.path_input_file)
    
