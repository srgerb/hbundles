#!/usr/bin/env python
from os import system,popen
import string
import argparse


#Takes in map of PDB file, separated into [row][column]
#Returns a list of chain lengths
#This method will work for any number of chains
def get_seg_lengths(pdb, start, end):
    count = 0 # current count of chain length
    chain = "A" # current chain name
    seg_length = 0 # List of lengths of chain in order
    for i in (0, len(input)):
        if input[i][4] != chain:
            chain = input[i][4]
            chain_lengths.append(count)
            count = 0
        count += 1
    return chain_lengths


#Swaps the ordering of the chains.
#Goes through the entire list for every chain. There is likely a more elegant way to do this. Mainly by going though the list once, counting the length of all the chains, and doing math to only iterate over the section you're pulling out, but this works for my purposes so far, and seems roughly the same time complexity for 2 chains
def swap_chains(chains):
    swapped = []
    for i in range(0, len(chains)):
        if len(chains[i]) >= 4:
            if chains[i][4] == "B":
                swapped.append(chains[i])
    for i in range(0, len(chains)):
        if len(chains[i]) >= 4:
            if chains[i][4] == "A":
                swapped.append(chains[i])
    renumber(swapped)
    return swapped

#takes in a list with information from a pdb in [row][column] format and renumbers the atoms, residues, and chain number
#need to figure out a more generalizable way to set up chain numbers.
def renumber(pdb):
    residue_number = 0 #count of residue number, PDB 6th column
    atom_number = 1 #atom number, pdb 2nd column
    prev_chain_num = 0
    prev_chain = pdb[0][4]
    chains = ["A","B","C","D","E","F","G","H","I","J"]
    for i in range(0, len(pdb)):
        if pdb[i][0] == "ATOM": #only modify Atom lines
            pdb[i][1] = atom_number
            atom_number += 1
            if pdb[i][2] == "N": #increments res number when it reaches the next nitrogen
                residue_number += 1
            pdb[i][5] = residue_number
            if pdb[i][4] != prev_chain: 
                prev_chain_num += 1
		prev_chain = pdb[i][4]
            else:
                pdb[i][4] = chains[prev_chain_num]

#takes in a list with information from a pdb in [row][column] format to appending from, a list to append to, int of the starting residue, and int of the end residue that is being appended
def append_to_list(list_from, list_to, start_res, end_res, rename_chain = None):
    for i in range(0, len(list_from)):
        if list_from[i][0] == "ATOM": #only add Atom lines
            residue_number = list_from[i][5]
            if int(residue_number) >= start_res:
                if int(residue_number) <= end_res:
                    if rename_chain is not None:
                        list_from[i][4] = rename_chain
                    list_to.append(list_from[i])
                else:
                    break

#takes 2 pdbs and writes pdb1 from residue 1 to first outer_seg_cutpoints, then pdb2 from first inner_seg_cutpoints to second inner_seg_cutpoints, then pdb1 from second outer_seg_cutpoints to the end of pdb1.
#this will only actually chimerize your protein if the proteins were aligned at the cutpoints
#need to set chains before doing this
def chimerize(pdb1, pdb2, outer_seg_cutpoints, inner_seg_cutpoints):
    chimera = []
    append_to_list(pdb1, chimera, 0, outer_segment_cutpoints[0], "A")
    append_to_list(pdb2, chimera, inner_seg_cutpoints[0], inner_seg_cutpoints[1],)
    append_to_list(pdb1, chimera, outer_segment_cutpoints[1], len(pdb1), "B")
    renumber(chimera)
    return chimera


#################################
###            MAIN           ###
#################################
#Purpose is to take in two PDB's and splice them together
#currently only works for chains that are perfectly aligned at splice sites

parser = argparse.ArgumentParser(description="Takes in 2 PDBs and their splice sites Returns PDB spliced together oSeg1-iSeg-oSeg2")
parser.add_argument('-s', type=str, nargs=2, help="Takes in 2 PDB files. The first will provide the outer segments of the output file (residues 1-w and x-end) and the second will provide the center segment of the output (residues x-z)")
parser.add_argument('-oSeg', type=int, nargs=2, default=(21,51), help="Specify outer segments. Output will be [1->arg1]-iSeg-[arg2->end]")
parser.add_argument('-iSeg', type=int, nargs=2, default=(20,58), help="Specify inner segments. Output will be [1->arg1]-iSeg-[arg2->end]")
parser.add_argument('-o', type=str, nargs=1, default="chimera", help="Name of output file")
parser.add_argument('--oFlip', type=bool, default=False, help="NOT RECCOMENDED! Will re-order the chains of the first (outer) pdb so that A is B and B is A. This changes the numbering.")
parser.add_argument('--iFlip', type=bool, default=False, help="NOT RECCOMENDED! Will re-order the chains of the second(inner) pdb so that A is B and B is A. This changes the numbering.")


args = parser.parse_args()
input_files=args.s #names of PDBs to splice together, in a list
outer_segment_cutpoints = args.oSeg
inner_segment_cutponits = args.iSeg
tag = args.o[0] #name of the output pdb

pdb_outer = map(string.split,open(input_files[0],'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]
pdb_inner = map(string.split,open(input_files[1],'r').readlines()) #open(input_file,'r') opens the input PDB in read mode, .readlines() function reads each line of that PDB, and map() applies string.split function to each of the lines, so input should be a 2D list of type [line][column]

#If oFlip or iFlip tags were used, swap chains A and B in that file
if args.oFlip == True:
    pdb_outer = swap_chains(pdb_outer)
if args.iFlip == True:
    pdb_inner = swap_chains(pdb_inner)

output = []
output = chimerize(pdb_outer, pdb_inner, outer_segment_cutpoints, inner_segment_cutponits)

#################################
###          OUTPUT           ###
#################################
output_file_name = tag
full_pdb = open('%s'%(output_file_name),'w')
for line in range(0, len(output)):
    if output[line][0] == 'ATOM':
        atom_name = output[line][2]
        if atom_name[0].isdigit():
            full_pdb.write('ATOM {atom_num:>6} {atom_type:4} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8} {occupancy} {temp_factor}\n'.format(atom_num = output[line][1],atom_type = output[line][2],residue = output[line][3], chain = output[line][4],res_num = output[line][5],coord_x = output[line][6],coord_y = output[line][7],coord_z = output[line][8],occupancy = output[line][9],temp_factor = output[line][10]))
        else:
            full_pdb.write('ATOM {atom_num:>6}  {atom_type:3} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8} {occupancy} {temp_factor}\n'.format(atom_num = output[line][1],atom_type = output[line][2],residue = output[line][3], chain = output[line][4],res_num = output[line][5],coord_x = output[line][6],coord_y = output[line][7],coord_z = output[line][8],occupancy = output[line][9],temp_factor = output[line][10]))
full_pdb.write('END')
