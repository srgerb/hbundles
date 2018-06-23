#!/usr/bin/env python
from os import system,popen
import string
import argparse

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        #if input_file[i] == ".":
        #    break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    return tag

def write_all(input_file,first_seg,second_seg,tag):
    pdb = map(string.split,open(input_file,'r').readlines())
    count1 = first_seg[0]
    while count1 <= first_seg[1]:
        count2  = second_seg[0]
        while count2 <= second_seg[1]:
            filename = tag+'_'+str(count1)+'_'+str(count2)+'.pdb'
            print (filename)
            write_pdb(pdb,count1,count2,filename)
            count2 +=1
        count1 += 1

def write_pdb(output,cut1,cut2,output_file_name):
    full_pdb = open('%s'%(output_file_name),'w')
    nterm = True
    for line in range(0, len(output)):
        if len(output[line]) >= 5:
            if output[line][0] == 'ATOM':
                if int(output[line][5]) < cut1 or int(output[line][5]) > cut2:
                    atom_name = output[line][2]
                    if atom_name[0].isdigit():
                        full_pdb.write('ATOM {atom_num:>6} {atom_type:4} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}\n'.format(atom_num = output[line][1],atom_type = output[line][2],residue = output[line][3], chain = output[line][4],res_num = output[line][5],coord_x = output[line][6],coord_y = output[line][7],coord_z = output[line][8]))
                    else:
                        full_pdb.write('ATOM {atom_num:>6}  {atom_type:3} {residue} {chain} {res_num:>3}    {coord_x:>8}{coord_y:>8}{coord_z:>8}\n'.format(atom_num = output[line][1],atom_type = output[line][2],residue = output[line][3], chain = output[line][4],res_num = output[line][5],coord_x = output[line][6],coord_y = output[line][7],coord_z = output[line][8]))
                if nterm == True and int(output[line][5]) > cut1: 
                    full_pdb.write("TER\n")
                    nterm = False
    full_pdb.write('END')

#takes in a list with information from a pdb in [row][column] format to appending from, a list to append to, int of the starting residue, and int of the end residue that is being appended
#I literally could just do this in the output area, which would make this marginally faster, but I'm modifying an old script and I'm lazy
def append_to_list(list_from, list_to, cutpoints, rename_chain = None):
    nterm = True
    for i in range(0, len(list_from)):
        if list_from[i][0] == "ATOM": #only add Atom lines
            residue_number = list_from[i][5]
        if int(residue_number) <= cutpoints[0]:
            if rename_chain is not None:
                list_from[i][4] = rename_chain
            list_to.append(list_from[i])
        elif (nterm):# this is just to insert TER in the cut site. this feels hacky and awful.
            ter = []
            ter.append("TER")
            list_to.append(ter)
            nterm = False
        if int(residue_number) >= cutpoints[1]:
            if rename_chain is not None:
                list_from[i][4] = rename_chain
            list_to.append(list_from[i])


#################################
###            MAIN           ###
#################################
#Purpose is to take in two PDB's and splice them together
#currently only works for chains that are perfectly aligned at splice sites

parser = argparse.ArgumentParser(description="Takes in PDBs two cut sites Returns cut PDB with TER listed in cut site")
parser.add_argument('PDB', type=str, nargs=1, help="PDB file to be spliced")
#parser.add_argument('cut', type=int, nargs=2, help="cut after first number and before second number")
parser.add_argument('-o', type=str, nargs=1, default=["cut"], help="Name of output file")

args = parser.parse_args()
tag = args.o[0] + '_' + clean_file_name(args.PDB[0]) #name of the output pdb

first_cut = clean_file_name(args.PDB[0]).split('_')
first_seg = (36,36)
cut_site = int(first_cut[7])
second_seg = (cut_site - 4, cut_site - 4)
print (second_seg)
#first_seg = (args.cut[0],args.cut[0]) #originally designed to cut in a range, between first and second number and print all of them
#second_seg = (args.cut[1],args.cut[1]) #I now just cuts at two points

write_all(args.PDB[0],first_seg,second_seg,tag)
