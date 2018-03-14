#!/usr/bin/env python
from os import system,popen
import string
import argparse

parser = argparse.ArgumentParser(description="Takes in a pdb with a pssm matrix in it and outputs that matrix.") 
parser.add_argument('PDB', type=str, nargs=1, help="PDB file to extract matrix from")
parser.add_argument('-o', type=str, nargs=1, default=['pssm'], help="Name of output file, default is 'pssm'")

args = parser.parse_args()
tag = args.o[0]

pdb_check = map(string.split,open(args.PDB[0],'r').readlines())
pdb = open(args.PDB[0],'r').readlines()
pssm = False
output = open('%s'%(tag),'w')
for i in range(0, len(pdb)):
    if (pdb_check[i] != None):
        if (pdb_check[i][0] == "#BEGIN_SEGMENT_SEQUENCE_PROFILE_PSSM"):
            pssm = True
        if pssm:
            output.write(pdb[i])
        if (pdb_check[i][0] == "#END_SEGMENT_SEQUENCE_PROFILE_PSSM"):
            pssm = False
            break
