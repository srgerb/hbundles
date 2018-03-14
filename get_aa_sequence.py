#!/usr/bin/env python
from os import system,popen
import string
import argparse

#takes a pdb file, and the remodel range, and fills nterm and cterm with the residue numbers and amino acid.
#nterm will be from start of pdb to remodel range i.e (1,S)..(range[0],N)
#cterm will be from end of remodel range to end of pdb i.e (range[1],F)..(120,G) if structure is 120 aa long.
def get_aa_sequence_from_pdb(pdb_file):
    pdb = map(string.split,open(pdb_file,'r').readlines())
    amino_acids = {"ALA":"A","GLY":"G","ILE":"I","LEU":"L","PRO":"P","VAL":"V","PHE":"F","TRP":"W","TYR":"Y","ASP":"D","GLU":"E","ARG":"R","HIS":"H","LYS":"K","SER":"S","THR":"T","CYS":"C","MET":"M","ASN":"N","GLN":"Q"}
    resnum = 0
    aa_sequence = []
    for i in range(0, len(pdb)):
        if (pdb[i][0] == "ATOM"):
            if (int(pdb[i][5]) != resnum):
                    resnum = int(pdb[i][5])
                    aa_sequence.append(amino_acids[pdb[i][3]])
    return aa_sequence

parser = argparse.ArgumentParser(description="Takes in a pdb and outputs a fasta file.") 
parser.add_argument('PDB', type=str, nargs=1, help="PDB file to be made into blueprint")
parser.add_argument('-o','--out', type=str, nargs=1, default="out", help="Name of output file")

args = parser.parse_args()
input_file=args.PDB[0] #name of PDB to extract the sequence from
tag = args.out[0] #name of the output pdb

output = get_aa_sequence_from_pdb(input_file)

output_file_name = tag + ".fasta"
fasta = open('%s'%(output_file_name),'w')
fasta.write(">" + tag + "\n")
for line in range(0, len(output)):
    fasta.write(output[i])
