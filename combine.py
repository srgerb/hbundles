#!/usr/bin/env python
from os import system,popen
import string
import argparse

#takes in a list with information from a pdb in [row][column] format and renumbers the atoms, residues, and chain number
#need to figure out a more generalizable way to set up chain numbers.
def renumber(pdb):
    residue_number = 0 #count of residue number, PDB 6th column
    atom_number = 1 #atom number, pdb 2nd column
    for i in range(1, len(pdb)):
        if pdb[i][0] == "ATOM" or pdb[i][0] == "HETATM": #only modify Atom lines
            pdb[i][1] = atom_number
            atom_number += 1
            if pdb[i][2] == "N": #increments res number when it reaches the next nitrogen
                residue_number += 1
            pdb[i][5] = residue_number

parser = argparse.ArgumentParser(description="Takes in 2 PDB files to be combined")
parser.add_argument('PDB', type=str, nargs=2, help=".pdb file to be renumbered")
parser.add_argument('-o', type=str, nargs=1, default=["out.pdb"], help="Name of output file")
args = parser.parse_args()
tag = args.o[0] #name of the output pdb

pdb = map(string.split,open(args.PDB[0],'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]
pdb += map(string.split,open(args.PDB[1],'r').readlines())
renumber(pdb)

connections = open(args.PDB[0],'r').readlines()
output = pdb
#################################
###          OUTPUT           ###
#################################
output_file_name = tag
full_pdb = open('%s'%(output_file_name),'w')
full_pdb.write(connections[0])
for line in range(1, len(output)):
    if output[line][0] == 'ATOM':
	atom_name = output[line][2]
        if atom_name[0].isdigit():
            full_pdb.write('ATOM {atom_num:>6} {atom_type:4} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}\n'.format(atom_num = output[line][1],atom_type = output[line][2],residue = output[line][3], chain = output[line][4],res_num = output[line][5],coord_x = output[line][6],coord_y = output[line][7],coord_z = output[line][8]))
        else:
            full_pdb.write('ATOM {atom_num:>6}  {atom_type:3} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}\n'.format(atom_num = output[line][1],atom_type = output[line][2],residue = output[line][3], chain = output[line][4],res_num = output[line][5],coord_x = output[line][6],coord_y = output[line][7],coord_z = output[line][8]))
