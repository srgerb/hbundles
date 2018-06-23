#!/usr/bin/env python
from os import system,popen
import string
import argparse

#cleans filename so "folder/input.pdb" will become "input"
def clean_file(input_file):
    tag = ""
    for i in range(0, len(input_file) - 4):# removes .pdb from end of file name
#        if input_file[i] == ".":
#            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    return tag

parser = argparse.ArgumentParser(description="Takes in 2 PDB to be renumbered call as 'python filt_chain.py $pdb $chain_letters' ")
parser.add_argument('PDB', type=str, nargs=1, help=".pdb file to be filtered")
parser.add_argument('chain', type=str, nargs='+',  help="Chains to be filtered out")
parser.add_argument('-o', type=str, nargs=1, help= "prefix to output file")
args = parser.parse_args()
input_file=args.PDB[0]
chain = args.chain #name of the output pdb

pdb = map(string.split,open(input_file,'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]

output = open(input_file,'r').readlines()
tag = clean_file(input_file)

for i in chain:
    tag += '_' + i
if args.o != None:
    output_file_name = args.o[0] + tag + ".pdb" #adds the filtered chain name to filename
else:
    output_file_name = tag + ".pdb"
pdb_chain = open('%s'%(output_file_name),'w')
for line in range(0, len(output)):
    if pdb[line] != []:
        if pdb[line][0] == 'TER':
            pdb_chain.write("TER\n")
        if pdb[line][0] == 'ATOM' and pdb[line][4] in chain:
            pdb_chain.write(output[line])
