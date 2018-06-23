from collections import namedtuple
import math
import statistics
import os

import pyrosetta

import argparse


parser = argparse.ArgumentParser(description="Takes in 2 PDB to be renumbered")
parser.add_argument('PDB', type=str, nargs=1, help=".pdb file to be filtered")
parser.add_argument('-o', type=str, nargs=1, help= "prefix to output file")
args = parser.parse_args()
pdb_file=args.PDB[0]

pyrosetta.init()
symmetry_mover = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover('/home/srgerb/hBundles/SymFiles/C4_Z.sym')
pose = pyrosetta.pose_from_file(pdb_file)
symmetry_mover.apply(pose)
tag = ""
for i in range(0, len(pdb_file) - 4):# removes .pdb from end of file name
            tag += pdb_file[i]
            if pdb_file[i] == "/": #removes path before file name
                tag = ""

pose.dump_pdb(tag + "_symm.pdb")
