#!/usr/bin/env python
from os import system,popen
import string
import argparse


parser = argparse.ArgumentParser(description="Takes in a score file and removes extra line")
parser.add_argument('PDB', type=str, nargs=1, help=".isc file to be fixed")
parser.add_argument('-o', type=str, nargs=1, default=["fixed.sc"], help="Name of output file")
args = parser.parse_args()
input_file=args.PDB[0]
tag = args.o[0] #name of the output pdb
columns = 68 #the number of columns your score file should have 
mismatch = 4 #the column that should be deleted

output = map(string.split,open(input_file,'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]

#################################
###          OUTPUT           ###
#################################
output_file_name = tag
full_pdb = open('%s'%(output_file_name),'w')
for line in range(0, len(output)):
    for column in range(0, len(output[line])):
        if (len(output[line]) != columns and column != mismatch):
            full_pdb.write(output[line][column] + "\t")
    full_pdb.write("\n")
