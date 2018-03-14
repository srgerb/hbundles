#!/usr/bin/env python
import os
from os import system,popen
import string
import glob


def check_dir(dirs, resnum):
    for i in dirs:
        if int(resnum) == int(i):
            return True
    return False
#Takes in map of PDB file, separated into [row][column]
#this does not scale well
def get_res_num(pdb_file, res_nums, output_file):
    pdb_map= map(string.split,open(pdb_file,'r').readlines())
    res_num = 0
    for line in pdb_map:
         if line[0] == 'ATOM' and int(line[5]) > int(res_num):
             res_num = line[5]
         elif line == None or line[0] == 'TER':
             break
    if check_dir(res_nums, res_num) == False:
        output_file.write("mkdir " + res_num + "\n")
        print("mkdir " + res_num)
        res_nums.append(res_num)
    output_file.write("mv " + pdb_file +" "+ res_num + "/ \n")
#    print("mv " + pdb_file +" "+ res_num)


#################################
###            MAIN           ###
#################################
res_nums = ['0']
output_file = open('%s'%('move.sh'),'a')

for i, pdb_file in enumerate(glob.glob('*.pdb')):
    get_res_num(pdb_file, res_nums, output_file)

