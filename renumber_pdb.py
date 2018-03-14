#!/usr/bin/env python
from os import system,popen
import string
import argparse

#takes in a list with information from a pdb in [row][column] format and renumbers the atoms, residues, and chain number
#need to figure out a more generalizable way to set up chain numbers.
def renumber(pdb):
    residue_number = 0 #count of residue number, PDB 6th column
    atom_number = 1 #atom number, pdb 2nd column
    prev_chain_num = pdb[0][4]
    for i in range(0, len(pdb)):
        if len(pdb[i]) > 4:
            if pdb[i][0] == "ATOM": #only modify Atom lines
                pdb[i][1] = atom_number
                atom_number += 1
                if pdb[i][2] == "N": #increments res number when it reaches the next nitrogen
                    residue_number += 1
                pdb[i][5] = residue_number
#                if pdb[i][4] == prev_chain_num: #hacky and awful, only works for two chain swap.
#                    pdb[i][4] = "A"
#                else:
#                    pdb[i][4] = "B"

parser = argparse.ArgumentParser(description="Takes in PDB to be renumbered")
parser.add_argument('PDB', type=str, nargs=1, help=".pdb file to be renumbered")
parser.add_argument('-o', type=str, nargs=1, default=["renumbered.pdb"], help="Name of output file")
args = parser.parse_args()
input_file=args.PDB[0]
tag = args.o[0] #name of the output pdb

pdb = map(string.split,open(input_file,'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]
renumber(pdb)
print ("renumbered")
connections = open(input_file,'r').readlines()
output = pdb
#################################
###          OUTPUT           ###
#################################
output_file_name = tag
full_pdb = open('%s'%(output_file_name),'w')
print ("Start output")
for line in output:
    if len(line) > 4:
        if line[0] == 'ATOM':
	    atom_name = line[2]
            if atom_name[0].isdigit():
                full_pdb.write('ATOM {atom_num:>6} {atom_type:4} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}  {occupancy}  {temp_factor}\n'.format(atom_num = line[1],atom_type = line[2],residue = line[3], chain = line[4],res_num = line[5],coord_x = line[6],coord_y = line[7],coord_z = line[8],occupancy = line[9],temp_factor = line[10]))
            else:
                full_pdb.write('ATOM {atom_num:>6}  {atom_type:3} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}  {occupancy}  {temp_factor}\n'.format(atom_num = line[1],atom_type = line[2],residue = line[3], chain = line[4],res_num = line[5],coord_x = line[6],coord_y = line[7],coord_z = line[8],occupancy = line[9],temp_factor = line[10],))
full_pdb.write('TER\n')
print("Print Connections")
for i in range(0, len(output)):
    if len(output[i]) > 1:
        if output[i][0] == 'CONECT':
            full_pdb.write(connections[i])
