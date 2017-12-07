#!/usr/bin/env python
from __future__ import print_function
import sys, os
import argparse
from pytriqs.applications.dft.sumk_dft import *
from pytriqs.applications.dcore.typed_parser import TypedParser
from pytriqs.applications.dcore.dmft_core import DMFTCoreSolver, create_parser

parser = argparse.ArgumentParser(\
        prog='dcore.py',\
        description='.',\
        epilog='end',\
        usage = '$ dcore input.ini',\
        add_help= True)

parser.add_argument('path_input_file', \
                    action = 'store',\
                    default= None,    \
                    type=str, \
                    help = "input file name.")

args=parser.parse_args()
if(os.path.isfile(args.path_input_file) is False):
    print("Input file is not exist.")
    sys.exit()

# Set Default value
parser = create_parser()

#
# Parse keywords and store
#
parser.read(args.path_input_file)
params = parser.as_dict()

solver = DMFTCoreSolver(params["model"]["seedname"], params)

solver.solve(max_step=params["control"]["max_step"], output_file=params["model"]["seedname"]+'.out.h5', output_group='dmft_out')
