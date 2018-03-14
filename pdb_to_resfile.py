#!/usr/bin/env python
from os import system,popen
import string
import argparse

def get_aa_sequence_from_pdb(pdb_file):
    pdb = map(string.split,open(pdb_file,'r').readlines())
    amino_acids = {"ALA":"A","GLY":"G","ILE":"I","LEU":"L","PRO":"P","VAL":"V","PHE":"F","TRP":"W","TYR":"Y","ASP":"D","GLU":"E","ARG":"R","HIS":"H","LYS":"K","SER":"S","THR":"T","CYS":"C","MET":"M","ASN":"N","GLN":"Q","HID":"H","HIE":"H"}
    resnum = 0
    res = []
    aa_sequence = []
    chain = []
    for i in range(0, len(pdb)):
        if (len(pdb[i]) >4):
            if (pdb[i][0] == "ATOM"):
                if (int(pdb[i][5]) != resnum):
                    resnum = int(pdb[i][5])
                    res.append(resnum)
                    aa_sequence.append(amino_acids[pdb[i][3]])
                    chain.append(pdb[i][4])
    sequence = []
    sequence.append(res)
    sequence.append(aa_sequence)
    sequence.append(chain)
    return sequence

parser = argparse.ArgumentParser(description="Takes in a pdb and outputs a resfile with same amino acid sequence.") 
parser.add_argument('PDB', type=str, nargs=1, help="PDB file to be made into blueprint")
parser.add_argument('-o','--out', type=str, nargs=1, default=["resfile"], help="Name of output file")

args = parser.parse_args()
input_file=args.PDB[0] #name of PDB to extract the sequence from
tag = args.out[0] #name of the output pdb

output = get_aa_sequence_from_pdb(input_file)

output_file_name = tag
resfile = open('%s'%(output_file_name),'w')
resfile.write("NATAA\nstart\n\n")
for i in range(0, len(output[0])):
    resfile.write("{resnumber} {chain} PIKAA {res_name}\n".format(resnumber = output[0][i],chain = output[2][i], res_name = output[1][i]))
