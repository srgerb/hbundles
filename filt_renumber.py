#!/usr/bin/env python
from os import system,popen
import string
import argparse

#takes in a list with information from a pdb in [row][column] format and renumbers the atoms, residues, and chain number
#only works for up to 26 chains
def renumber(pdb, filt_chain):
    chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    residue_number = 0 #count of residue number, PDB 6th column
    atom_number = 1 #atom number, pdb 2nd column
    prev_chain = pdb[0][4] #this could cause errors if the first line is not ATOM
    chain_num = 0 #count of chain number, 0 ordered
    renumbered_pdb = []
    for i in range(0, len(pdb)):
        if len(pdb[i]) > 4:
            if pdb[i][4] not in filt_chain:
                break
            if pdb[i][0] == "ATOM": #only modify Atom lines
                pdb[i][1] = atom_number
                atom_number += 1
                if pdb[i][2] == "N": #increments res number when it reaches the next nitrogen
                    residue_number += 1
                    if i != 0:
                        if pdb[i][4] != prev_chain or pdb[i-1][0] == "TER":
                            chain_num += 1
                            prev_chain = pdb[i][4]
                pdb[i][4] = chains[chain_num]        
                pdb[i][5] = residue_number
        renumbered_pdb.append(pdb[i])
    return renumbered_pdb

#cleans filename so "folder/input.pdb" will become "input"
def clean_output_filename(input_file, prefix, chain):
    tag = ""
    for i in range(0, len(input_file) - 4):# removes .pdb from end of file name
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    chain_names = ""
    for i in chain:
        chain_names += i
    if prefix != None:
        output_file_name = prefix[0] + tag + "_" + chain_names + ".pdb" #adds the filtered chain name to filename
    else:
        output_file_name = tag + "_" + chain_names + ".pdb"
    return output_file_name

def write_file(renumbered_pdb_map, output_file_name):
    filtered_pdb = open('%s'%(output_file_name),'w')
    #print ("Start output")
    for line in renumbered_pdb_map:
        if len(line) > 4:
            if line[0] == 'ATOM':
                atom_name = line[2]
                if atom_name[0].isdigit():
                    filtered_pdb.write('ATOM {atom_num:>6} {atom_type:4} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}  {occupancy}  {temp_factor}\n'.format(atom_num = line[1],atom_type = line[2],residue = line[3], chain = line[4],res_num = line[5],coord_x = line[6],coord_y = line[7],coord_z = line[8],occupancy = line[9],temp_factor = line[10]))
                else:
                    filtered_pdb.write('ATOM {atom_num:>6}  {atom_type:3} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}  {occupancy}  {temp_factor}\n'.format(atom_num = line[1],atom_type = line[2],residue = line[3], chain = line[4],res_num = line[5],coord_x = line[6],coord_y = line[7],coord_z = line[8],occupancy = line[9],temp_factor = line[10],))
    filtered_pdb.write('TER\n')

parser = argparse.ArgumentParser(description="Takes in 2 PDB to be renumbered")
parser.add_argument('PDB', type=str, nargs=1, help=".pdb file to be filtered")
parser.add_argument('chain', type=str, nargs='+',  help="Chains to be filtered out")
parser.add_argument('-o', type=str, nargs=1, help= "prefix to output file")
args = parser.parse_args()
input_file=args.PDB[0]
chain = args.chain #name of the chains

pdb = map(string.split,open(input_file,'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]
renumbered_pdb = renumber(pdb, chain)
output_file_name = clean_output_filename(input_file, args.o, chain)
write_file(renumbered_pdb, output_file_name)
