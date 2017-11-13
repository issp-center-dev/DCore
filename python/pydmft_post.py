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

#from .typed_parser import TypedParser
from typed_parser import TypedParser


def __print_paramter(p, param_name):
    print(param_name + " = " + str(p[param_name]))

def __generate_wannier90_model(params, l, norb, equiv, n_k, kvec):
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
    n_k : integer
        Number of k points
    kvec : float array
        k-points where A(k,w) is computed

    Returns
    -------
    hopping : complex
        k-dependent one-body Hamiltonian
    n_orbitals : integer
        Number of orbitals at each k. It does not depend on k
    proj_mat : complex
        Projection onto each correlated orbitals
    """
    ncor = params["ncor"]
    n_spin = 1

    print("               ncor = ", ncor)
    for i in range(ncor):
        if equiv[i] == -1: equiv[i] = i
        print("     l[{0}], norb[{0}], equiv[{0}] = {1}, {2}, {3}".format(i,l[i],norb[i],equiv[i]))
    #
    # Read hopping in the real space from the Wannier90 output
    #
    w90c = Wannier90Converter(seedname=params["seedname"])
    nr, rvec, rdeg, nwan, hamr = w90c.read_wannier90hr(params["seedname"]+"_hr.dat")
    #
    # Fourier transformation of the one-body Hamiltonian
    #
    n_orbitals = numpy.ones([n_k, n_spin], numpy.int) * nwan
    hopping = numpy.zeros([n_k, n_spin, numpy.max(n_orbitals), numpy.max(n_orbitals)], numpy.complex_)
    for ik in range(n_k):
        for ir in range(nr):
            rdotk = 2 * numpy.pi * numpy.dot(kvec[ik,:], rvec[ir,:])
            factor = (numpy.cos(rdotk) + 1j * numpy.sin(rdotk)) / float(rdeg[ir])
            hopping[ik,0,:, :] += factor * hamr[ir][:, :]
    #
    # proj_mat is (norb*norb) identities at each correlation shell
    #
    proj_mat = numpy.zeros([n_k, n_spin, ncor, numpy.max(norb), numpy.max(n_orbitals)], numpy.complex_)
    iorb = 0
    for icor in range(ncor):
        proj_mat[:, :, icor, 0:norb[icor], iorb:iorb + norb[icor]] = numpy.identity(norb[icor], numpy.complex_)
        iorb += norb[icor]

    return hopping, n_orbitals, proj_mat

# FIXME: split this LONG function
def __generate_lattice_model(params, n_k, kvec):
    """
    Compute hopping etc. for A(k,w) of preset models

    Parameters
    ----------
    params : dictionary
        Input parameters
    n_k : integer
        Number of k points
    kvec : float array
        k-points where A(k,w) is computed

    Returns
    -------
    hopping : complex
        k-dependent one-body Hamiltonian
    n_orbitals : integer
        Number of orbitals at each k. It does not depend on k
    proj_mat : complex
        Projection onto each correlated orbitals
    """
    #
    # Construct model
    #
    if params["orbital_model"] == 'single':
        norb = 1
    elif params["orbital_model"] == 'eg':
        #FIXME: l=2 does not make sense. l=2 assumes norb=5 (full d-shell) in generating Coulomb tensor.
        #What is the proper way to generate Coulomb tensor for eg?
        norb = 2
    elif params["orbital_model"] == 't2g':
        norb = 3
    elif params["orbital_model"] == 'full-d':
        norb = 5
    else:
        print("Error ! Invalid lattice : ", params["orbital_model"])
        sys.exit()

    t = params["t"]
    tp = params["t'"]
    n_spin = 1
    #
    # Energy band
    #
    if params["lattice"] == 'bethe':
        #
        # For Bhete lattice, k-point has no meanings.
        #
        print("Skip")
    else:
        n_orbitals = numpy.ones([n_k, n_spin], numpy.int) * norb
        hopping = numpy.zeros([n_k, n_spin, norb, norb], numpy.complex_)

        for ik in range(n_k):
            if params["lattice"] == 'chain':
                ek = 2.0*t*numpy.cos(kvec[ik,0]) + 2*tp*numpy.cos(2.0*kvec[ik,0])
            elif params["lattice"] == 'square':
                ek = 2.0 * t * (numpy.cos(kvec[ik,0]) + numpy.cos(kvec[ik,1])) \
                   + 2.0 * tp * (numpy.cos(kvec[ik,0] + kvec[ik,1]) + numpy.cos(kvec[ik,0] - kvec[ik,1]))
            elif params["lattice"] == 'cubic':
                ek = 2 * t * (numpy.cos(kvec[ik,0]) + numpy.cos(kvec[ik,1]) + numpy.cos(kvec[ik,2])) \
                    + 2 * tp * (numpy.cos(kvec[ik,0] + kvec[ik,1]) + numpy.cos(kvec[ik,0] - kvec[ik,1]) \
                              + numpy.cos(kvec[ik,1] + kvec[ik,2]) + numpy.cos(kvec[ik,1] - kvec[ik,2]) \
                              + numpy.cos(kvec[ik,2] + kvec[ik,0]) + numpy.cos(kvec[ik,2] - kvec[ik,0]))

            for iorb in range(norb): hopping[ik, 0, iorb, iorb] = ek
    #
    # proj_mat is (norb*norb) identities at each correlation shell
    #
    proj_mat = numpy.zeros([n_k, n_spin, 1, norb, norb], numpy.complex_)
    proj_mat[:, :, 0, 0:norb, 0:norb] = numpy.identity(norb, numpy.complex_)

    return hopping, n_orbitals, proj_mat


def pydmft_post(filename):
    """
    Main routine for the post-processing tool

    Parameters
    ----------
    filename : string
        Input-file name
    """
    print("Reading {0} ...\n".format(filename))
    #
    # Construct a parser with default values
    #
    p = TypedParser();
    p.add_option("model", "t", float, 1.0, "some help message")
    p.add_option("model", "t'", float, 0.0, "some help message")
    p.add_option("model", "orbital_model", str, "single", "some help message")
    p.add_option("model", "nk", int, 8, "some help message")
    p.add_option("model", "ncor", int, 1, "some help message")
    p.add_option("model", "lattice", str, "chain", "some help message")
    p.add_option("model", "seedname", str, "pydmft", "some help message")
    p.add_option("model", "cshell", str, "[]", "some help message")
    p.add_option("model", "nnode", int, 2, "some help message")
    p.add_option("model", "knode", str, "[]", "some help message")
    p.add_option("model", "bvec", str, "[]", "some help message")
    p.read(filename)
    p_model = p.as_dict()["model"]
    #
    # Information of correlation shells. It is used only in conjunction to Wannier90.
    # cshell=(l, norb, equiv) or (l, norb)
    #
    cshell_list=re.findall(r'\(\s*\d+,\s*\d+,*\s*\d*\)', p_model["cshell"])
    l = [0]*p_model['ncor']
    norb = [1]*p_model['ncor']
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
    #
    # Nodes for k-point path
    # knode=(label,k0, k1, k2) in the fractional coordinate
    #
    knode_list = re.findall(r'\(\w\d?,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)', p_model["knode"])
    knode = numpy.zeros((p_model['nnode'], 3), numpy.float_)
    klabel =['G'] * p_model['nnode']
    try:
        for i, _list in enumerate(knode_list):
            _knode = filter(lambda w: len(w) > 0, re.split(r'[\(\s*\,\s*\,\s*,\s*\)]', _list))
            klabel[i] = _knode[0]
            for j in range(3): knode[i,j] = float(_knode[j+1])
    except:
        raise RuntimeError("Error ! Format of knode is wrong.")
    #
    # Reciprocal lattice vectors
    # bvec=[(b0x, b0y, k0z),(b1x, b1y, k1z),(b2x, b2y, k2z)]
    #
    bvec_list = re.findall(r'\(\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*,\s*-?\s*\d+\.?\d*\)', p_model["bvec"])
    bvec = numpy.zeros((3, 3), numpy.float_)
    try:
        for i, _list in enumerate(bvec_list):
            _bvec = filter(lambda w: len(w) > 0, re.split(r'[\(\s*\,\s*,\s*\)]', _list))
            for j in range(3): bvec[i,j] = float(_bvec[j])
    except:
        raise RuntimeError("Error ! Format of bvec is wrong.")
    #
    # Summary of input parameters
    #
    print("Parameter summary")
    for k,v in p_model.items():
        print(k, " = ", str(v))
    #
    # Construct parameters for the A(k,w)
    #
    nk_line = p_model["nk"]
    n_k = (p_model["nnode"] - 1)*nk_line + 1
    print(" Total number of k =", str(n_k))
    kvec = numpy.zeros((n_k, 3), numpy.float_)
    ikk = 0
    for inode in range(p_model["nnode"] - 1):
        for ik in range(nk_line + 1):
            if inode != 0 and ik == 0: continue
            for i in range(3):
                kvec[ikk,i] = float((nk_line - ik)) * knode[inode,i] + float(ik) * knode[inode + 1,i]
                #kvec[ikk,i] = 2.0 * numpy.pi * kvec[ikk,i] / float(nk_line)
                kvec[ikk, i] = kvec[ikk, i] / float(nk_line)
            ikk += 1

    if p_model["lattice"] == 'wannier90':
        hopping, n_orbitals, proj_mat = __generate_wannier90_model(p_model, l, norb, equiv, n_k, kvec)
    else:
        hopping, n_orbitals, proj_mat = __generate_lattice_model(p_model, n_k, kvec)
    #
    # Output them into seedname.h5
    #
    f = HDFArchive(p_model["seedname"]+'.h5','a')
    if not ("dft_bands_input" in f):
        f.create_group("dft_bands_input")
    f["dft_bands_input"]["hopping"] = hopping
    f["dft_bands_input"]["n_k"] = n_k
    f["dft_bands_input"]["n_orbitals"] = n_orbitals
    f["dft_bands_input"]["proj_mat"] = proj_mat
    del  f
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
    
