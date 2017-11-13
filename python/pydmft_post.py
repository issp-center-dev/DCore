#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import numpy
import math
import argparse
import re
from pytriqs.archive.hdf_archive import HDFArchive
from pytriqs.applications.dft.converters.wannier90_converter import Wannier90Converter
from pytriqs.applications.dft.converters.hk_converter import HkConverter

#from .typed_parser import TypedParser
from typed_parser import TypedParser


def __print_paramter(p, param_name):
    print(param_name + " = " + str(p[param_name]))

def __generate_wannier90_model(params, l, norb, equiv, f, knode):
    nk = params["nk"]
    nk0 = params["nk0"]
    nk1 = params["nk1"]
    nk2 = params["nk2"]
    ncor = params["ncor"]
    nelec = params["nelec"]

    if nk0 == 0: nk0 = nk
    if nk1 == 0: nk1 = nk
    if nk2 == 0: nk2 = nk
    print("                nk0 = ", nk0)
    print("                nk1 = ", nk1)
    print("                nk2 = ", nk2)
    print("               ncor = ", ncor)
    for i in range(ncor):
        if equiv[i] == -1: equiv[i] = i
        print("     l[{0}], norb[{0}], equiv[{0}] = {1}, {2}, {3}".format(i,l[i],norb[i],equiv[i]))

    print("0 {0} {1} {2}".format(nk0,nk1,nk2), file=f)
    print(nelec, file=f)
    print(ncor, file=f)
    for i in range(ncor):
        print("{0} {1} {2} {3} 0 0".format(i,equiv[i],l[i],norb[i]), file=f)

    #w90 = Wannier90Converter(seedname = seedname)
    nr, rvec, rdeg, nw, hamr = Wannier90Converter.read_wannier90hr(params["seedname"]+"_hr.dat")

    kvec = [0.0, 0.0, 0.0]
    twopi = 2 * numpy.pi
    h_of_k = [numpy.zeros((nw, nw), dtype=numpy.complex_)
              for ik in range((params["nnode"] - 1)*params["nk"]+1)]
    for inode in range(params["nnode"] - 1):
        for ik in range(nk + 1):
            if inode != 0 & ik == 0: pass
            for i in range(3):
                kvec[i] = ik * knode[i + 3 * inode] + (nk - ik) * knode[i + 3 * (inode + 1)]
                kvec[i] = 2.0 * numpy.pi * kvec[i] / float(nk)

                for ir in range(nr):
                    rdotk = twopi * numpy.dot(kvec, rvec[ir])
                    factor = (math.cos(rdotk) + 1j * math.sin(rdotk)) / \
                             float(rdeg[ir])
                    h_of_k[ik][:, :] += factor * hamr[ir][:, :]

# FIXME: split this LONG function
def __generate_lattice_model(params, l, norb, equiv, f, knode):
    __print_paramter(params, "t")
    __print_paramter(params, "tp")
    __print_paramter(params, "nk")
    __print_paramter(params, "orbital_model")
    weights_in_file = False
    if params["lattice"] == 'chain':
        nkBZ = params["nk"]
    elif params["lattice"] == 'square':
        nkBZ = params["nk"]**2
    elif params["lattice"] == 'cubic':
        nkBZ = params["nk"]**3
    elif params["lattice"] == 'bethe':
        nkBZ = params["nk"]
        weights_in_file = True
    else:
        print("Error ! Invalid lattice : ", params["lattice"])
        sys.exit()
    print(" Total number of k =", str(nkBZ))

    #
    # Model
    #
    if params["orbital_model"] == 'single':
        l[0] = 0
        norb[0] = 1
    elif params["orbital_model"] == 'eg':
        #FIXME: l=2 does not make sense. l=2 assumes norb=5 (full d-shell) in generating Coulomb tensor.
        #What is the proper way to generate Coulomb tensor for eg?
        l[0] = 2
        norb[0] = 2
    elif params["orbital_model"] == 't2g':
        l[0] = 1
        norb[0] = 3
    elif params["orbital_model"] == 'full-d':
        l[0] = 2
        norb[0] = 5
    else:
        print("Error ! Invalid lattice : ", params["orbital_model"])
        sys.exit()
    #
    # Write General-Hk formatted file
    #
    print(nkBZ, file=f)
    print(params["nelec"], file=f)
    print("1", file=f)
    print("0 0 {0} {1}".format(l[0], norb[0]), file=f)
    print("1", file=f)
    print("0 0 {0} {1} 0 0".format(l[0], norb[0]), file=f)
    print("1 {0}".format(norb[0]), file=f)

    t = params["t"]
    tp = params["tp"]
    nk = params["nk"]

    #
    # Energy band
    #
    kvec = [0.0, 0.0, 0.0]
    if params["lattice"] == 'bethe':
        #
        # If Bethe lattice, set k-weight manually to generate semi-circular DOS
        #
        print("Skip")
    else:
        for inode in range(params["nnode"] - 1):
            for ik in range(nk+1):
                if inode != 0 & ik == 0: pass
                for i in range(3):
                    kvec[i] = ik * knode[i+3*inode] + (nk-ik)*knode[i+3*(inode+1)]
                    kvec[i] = 2.0 * numpy.pi * kvec[i] / float(nk)

                    if params["lattice"] == 'chain':
                        ek = 2.0*t*numpy.cos(kvec[0]) + 2*tp*numpy.cos(2.0*kvec[0])
                    elif params["lattice"] == 'square':
                        ek = 2.0 * t * (numpy.cos(kvec[0]) + numpy.cos(kvec[1])) \
                           + 2.0 * tp * (numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1]))
                    elif params["lattice"] == 'cubic':
                        ek = 2 * t * (numpy.cos(kvec[0]) + numpy.cos(kvec[1]) + numpy.cos(kvec[2])) \
                           + 2 * tp * (numpy.cos(kvec[0] + kvec[1]) + numpy.cos(kvec[0] - kvec[1]) \
                                     + numpy.cos(kvec[1] + kvec[2]) + numpy.cos(kvec[1] - kvec[2]) \
                                     + numpy.cos(kvec[2] + kvec[0]) + numpy.cos(kvec[2] - kvec[0]))

                    for iorb in range(norb[0]):
                        for jorb in range(norb[0]):
                            if iorb == jorb:
                                print("{0}".format(ek), file=f) #Real part
                            else:
                                print("0.0", file=f) #Real part
                    for iorb in range(norb[0]*norb[0]): print("0.0", file=f) #Imaginary part
    return weights_in_file

def pydmft_post(filename):
    print("Reading {0} ...\n".format(filename))
    #
    # Construct a parser with default values
    #
    p = TypedParser();
    p.add_option("model", "t", float, 1.0, "some help message")
    p.add_option("model", "tp", float, 0.0, "some help message")
    p.add_option("model", "orbital_model", str, "single", "some help message")
    p.add_option("model", "nk", int, 8, "some help message")
    p.add_option("model", "ncor", int, 1, "some help message")
    p.add_option("model", "lattice", str, "chain", "some help message")
    p.add_option("model", "seedname", str, "pydmft", "some help message")
    p.add_option("model", "cshell", str, "[]", "some help message")
    p.add_option("model", "nnode", int, 1, "some help message")
    p.add_option("model", "knode", str, "[]", "some help message")
    p.add_option("model", "bvec", str, "[]", "some help message")
    p.read(filename)
    p_model = p.as_dict()["model"]

    #cshell=(l, norb, equiv) or (l, norb)
    cshell_list=re.findall(r'\(\s*\d+,\s*\d+,*\s*\d*\)', p_model["cshell"])
    print("debug4 {0}\n {1}".format(cshell_list, p_model["cshell"]))
    l = [0]*p_model['ncor']
    norb = [0]*p_model['ncor']
    equiv = [-1]*p_model['ncor']
    try:
        for  i, _list  in enumerate(cshell_list):
            _cshell = filter(lambda w: len(w) > 0, re.split(r'[\(\s*\,\s*,*\s*\)]', _list))
            l[i] = int(_cshell[0])
            norb[i] = int(_cshell[1])
            if len(_cshell)==3:
                equiv[i] = int(_cshell[2])
    except:
        raise RuntimeError("Error ! Format of cshell is wrong.")

    # knode=(label,k0, k1, k2)
    knode_list = re.findall(r'\(\w\d?,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)', p_model["knode"])
    print("debug3 {0}".format(knode_list))
    knode = [0.0] * 3 * p_model['nnode']
    klabel =['G'] * p_model['nnode']
    try:
        for i, _list in enumerate(knode_list):
            _knode = filter(lambda w: len(w) > 0, re.split(r'[\(\s*\,\s*\,\s*,\s*\)]', _list))
            klabel[i] = _knode[0]
            for j in range(3): knode[j+i*3] = float(_knode[j+1])
    except:
        raise RuntimeError("Error ! Format of knode is wrong.")

    for i in range(p_model['nnode']):
        print("debug2 {0} {1} {2} {3}".format(klabel[i], knode[0+i*3],knode[1+i*3],knode[2+i*3]))
    # bvec=[(b0x, b0y, k0z),(b1x, b1y, k1z),(b2x, b2y, k2z)]
    bvec_list = re.findall(r'\(\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)', p_model["bvec"])
    bvec = [0.0] * 9
    try:
        for i, _list in enumerate(bvec_list):
            _bvec = filter(lambda w: len(w) > 0, re.split(r'[\(\s*\,\s*,\s*\)]', _list))
            for j in range(3): bvec[j+i*3] = float(_bvec[j])
    except:
        raise RuntimeError("Error ! Format of bvec is wrong.")
    for i in range(3):
        print("debug1 {0} {1} {2}".format(bvec[0+i*3],bvec[1+i*3],bvec[2+i*3]))

    #
    # Summary of input parameters
    #
    print("Parameter summary")
    for k,v in p_model.items():
        print(k, " = ", str(v))

    #
    # Lattice
    #
    seedname = p_model["seedname"]
    if p_model["lattice"] == 'wannier90':
        with open(seedname+'.inp', 'w') as f:
            __generate_wannier90_model(p_model, l, norb, equiv, f, knode)
        # Convert General-Hk to SumDFT-HDF5 format
    #    Converter = Wannier90Converter(seedname = seedname)
    #    Converter.convert_dft_input()
    else:
        with open(seedname+'.inp', 'w') as f:
            weights_in_file = __generate_lattice_model(p_model, l, norb, equiv, f, knode)
        # Convert General-Hk to SumDFT-HDF5 format
    #    Converter = HkConverter(filename = seedname + ".inp", hdf_filename=seedname+".h5")
    #    Converter.convert_dft_input(weights_in_file=weights_in_file)

    #
    # Finish
    #
    print("Done")
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(\
        prog='pydmft_post.py',\
        description='pre script for pydmft.',\
        epilog='end',\
        usage = '$ pydmft_post input',\
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
    pydmft_post(args.path_input_file)
    
