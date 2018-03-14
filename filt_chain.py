#!/usr/bin/env python
from os import system,popen
import string
import argparse


parser = argparse.ArgumentParser(description="Takes in 2 PDB to be renumbered")
parser.add_argument('PDB', type=str, nargs=1, help=".pdb file to be filtered")
parser.add_argument('chain', type=str, nargs=1,  help="Chain to be filtered out")
parser.add_argument('-o', type=str, nargs=1, help= "prefix to output file")
args = parser.parse_args()
input_file=args.PDB[0]
chain = args.chain[0] #name of the output pdb

pdb = map(string.split,open(input_file,'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]

output = open(input_file,'r').readlines()
tag = ""

#cleans filename so "folder/input.pdb" will become "input"
for i in range(0, len(input_file) - 4):# removes .pdb from end of file name
    tag += input_file[i]
    if input_file[i] == "/": #removes path before file name
       tag = ""
if args.o != None:
    output_file_name = args.o[0] + tag + "_" + chain + ".pdb" #adds the filtered chain name to filename
else:
    output_file_name = tag + "_" + chain + ".pdb"
pdb_chain = open('%s'%(output_file_name),'w')
for line in range(0, len(output)):
    if len(pdb[line]) > 4:
        if pdb[line][0] == 'ATOM' and pdb[line][4] == chain:
            pdb_chain.write(output[line])
