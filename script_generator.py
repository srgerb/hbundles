#!/usr/bin/env python
from os import system,popen
import string
import argparse

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        if input_file[i] == ".":
            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    return tag

def get_HID_resnum(input_file, commands_file):
    pdb = map(string.split,open(input_file,'r').readlines())
    for line in pdb:
        if (line != None and line[0] == "ATOM"):
            if line[3] == "HID":
                return line[5]
            elif line[0] == "TER" or line[0] == "END":
                print("No HID on first chain")
                return "0"
    print("No HID in pose")
    return "0"

def get_HID_and_HBNet(input_file,HID_resnum,hbnet_residues, commands_file):
    pdb = map(string.split,open(input_file,'r').readlines())
    for line in pdb:
        if len(line) != 0:
            if line[0] == "ATOM":
                if line[3] == "HID":
                    HID_resnum = line[5]
                elif line[0] == "TER" or line[0] == "END":
                    print("No HID on first chain")
            if line[0] == "REMARK":
                hbnet_residues.append(line[2])
    hbnet_var = ""
    for i in range(0, len(hbnet_residues)-1):
        hbnet_var += hbnet_residues[i] + ","
    hbnet_var += hbnet_residues[len(hbnet_residues)-1]
    commands_file.write("/home/srgerb/Rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease @flags -parser:script_vars hbnet_res=" + hbnet_var + " binder=" + HID_resnum + " constraints=/home/srgerb/hBundles/pl061_backbone_library/6_second_shell_interactions/csts/" + HID_resnum + ".cst -in:file:s " + input_file)

def get_cst_file_name(input_file, commands_file):
    tag= ""
    for i in range(0, len(input_file) - 9):# removes .pdb from end of file name
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    cst_file= "/home/srgerb/hBundles/pl061_backbone_library/6e_second_shell_interactions/renamed_csts/" + tag + ".cst"
    commands_file.write("/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/main/database/ @flags  -parser:script_vars csts=" + cst_file + " -in:file:s " + input_file)

def print_basic(input_file, commands_file):
    commands_file.write("/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/main/database/ @flags -in:file:s " + input_file)

def print_length(input_file, commands_file):
    pdb_map= map(string.split,open(input_file,'r').readlines())
    res_num = 0
    for line in pdb_map:
         if line[0] == 'ATOM' and int(line[5]) > int(res_num):
             res_num = line[5]
         elif line == None or line[0] == 'TER':
             break
    commands_file.write("/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/main/database/ @flags -parser:script_vars length=" + res_num + " -in:file:s " + input_file + "\n\n")

def make_folder_and_print_length(input_file, commands_file):
    tag = clean_file_name(input_file)
    output_file = open('%s'%('make_files.sh'),'a')
    output_file.write("mkdir " + tag + "\n")
    commands_file.write("cd " + tag + " && ")
    print_length(input_file, commands_file)

parser = argparse.ArgumentParser(description="Takes in a file and generates command to run it")
parser.add_argument('PDB', type=str, nargs=1, help=".pdb file that was generated by symm matcher")
args = parser.parse_args()
input_file = args.PDB[0]

commands_file = open('%s'%('array_tasks.list'),'a')
make_folder_and_print_length(input_file, commands_file)

#if args.o != None:
#    output_file_name = args.o[0] + "array_tasks.list" #adds the filtered chain name to filename
#else:
#    output_file_name = "array_tasks.list"

#cstfile = open('%s'%(output_file_name),'w')
#for i in range(0, last_chain_number):
#    cstfile.write("\n".format(resnumber = HID_resnum, chain = chains[i], next_chain = chains[(i+1)%last_chain_number] ))
#    cstfile.write("\n")

