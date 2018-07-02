import glob             #inputfiles (maybe)
import argparse         # parse arguments from command line
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import savefig

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        if input_file[i] == ".":
            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
            tag = ""
    return tag

def get_dssp(pdb_file):
    import pyrosetta 
    pyrosetta.init('-out:file:output_virtual true')
    pose = pyrosetta.pose_from_file(pdb_file)
    struct_from_pose = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
    dssp_secstruct =struct_from_pose.get_dssp_secstruct()
    print(dssp_secstruct)
    return dssp_secstruct

#make a list of tuples with start and end of each loop
def get_loops(dssp_secstruct):
    in_loop = False
    structure = False
    start_position = 0
    end_position = 0
    loops = []
    for i, residue_secstruct in enumerate(dssp_secstruct):
        if residue_secstruct != 'L' and not in_loop:
            structure = True
            continue
        elif residue_secstruct == 'L' and not in_loop and structure:
            if i >= 2:
                start_position = i - 2 # include structured residues likely part of the beta bulge
            elif i == 1:
                start_position = i - 1
            in_loop = True
            structure = False
        elif residue_secstruct == 'L' and in_loop:
            end_position = i + 1 # include structured residues likely part of the beta bulge
        elif residue_secstruct != 'L' and in_loop:
            loops.append((start_position,end_position))
            in_loop = False;
            structure = True
    return loops

#make a single list of all the residues
def get_all_res(loops):
    all_res = {}
    for loopset in loops:
        for i in range(loopset[0], loopset[1] + 1):
            all_res[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
    return all_res


def compare(fasta_name,all_res, loops):
    fasta = open(fasta_name, "r")
    for line in fasta.readlines():
        if line != [] and line[0] != '>':
            for i,residue in enumerate(line):
                if i in all_res:
                    #print ("increment " + str(i) + " "+ residue)
                    if residue in all_res[i]:
                        all_res[i][residue] += 1
    for loopset in loops:
        all_res[loopset[1] + 2] = {"A":1000,"C":1000,"D":1000,"E":1000,"F":1000,"G":1000,"H":1000,"I":1000,"K":1000,"L":1000,"M":1000,"N":1000,"P":1000,"Q":1000,"R":1000,"S":1000,"T":1000,"V":1000,"W":1000,"Y":1000}


def make_dataframe(all_res,fasta_name):
    dataframe = pd.DataFrame(all_res)
    heatmap = sns.heatmap(dataframe, cmap="YlGnBu", xticklabels = 2)
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=90)
    fig = heatmap.get_figure()
    plot_name = clean_file_name(fasta_name) + ".png"
    fig.savefig(plot_name)
    return dataframe 
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Takes in any number of pdbs and searches through to find a symmetric metal site based on constraints in this file.")
    parser.add_argument('PDBS', type=str, nargs=1, help="pdbs to be searched through. Can be /path/to/files/*.pdb or /path/to/your_file.pdb or your.pdb etc")
    parser.add_argument('Fasta', type=str, nargs=1, help="fasta file of aligned proteins to compare")
    parser.add_argument('-o','--out', type=str, nargs=1, default=[""], help="prefix for the output file")
                    
    args = parser.parse_args()
    input_files=args.PDBS[0] #name of PDB to extract the sequence from
    tag = args.out[0] #name of the output pdb

                                    
    for i, pdb_file in enumerate(glob.glob(input_files)):
        dssp_string = get_dssp(pdb_file)
        # dssp_string = "HLLLLLLLLH"
        loops = get_loops(dssp_string)
        print (loops)
        all_res = get_all_res(loops)
        # print ("all_res[4]")
        # print (all_res[4])
        # print ("all_res[4][A]")
        # print (all_res[4]['A'])
        #print(all_res.values())
        compare(args.Fasta[0],all_res, loops)
        # print ("all_res[4]")
        # print (all_res[4])
        # print ("all_res[4][A]")
        # print (all_res[4]['A'])
        # print ("all_res[5]")
        # print (all_res[5])
        # print ("all_res[5][A]")
        # print (all_res[5]['A'])
        make_dataframe(all_res, args.Fasta[0])
        # make_heatmap(all_res)
