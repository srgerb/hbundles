import glob             #inputfiles (maybe)
import argparse         # parse arguments from command line
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import savefig
import re

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        if input_file[i] == ".":
            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
            tag = ""
    return tag

#make a list of tuples with start and end of each loop 
#dssp_secstruct is a string of the secondary structure of the protein
#ss_desired is a string letter for the secondary structure desired('E' for strand, 'H' for helix, 'L' for loop)
def get_structure(dssp_secstruct, ss_desired, min_length = 0, max_length = -1):
    in_loop = False
    structure = False
    start_position = 0
    end_position = 0
    loops = []
    for i, residue_secstruct in enumerate(dssp_secstruct):
        if residue_secstruct != ss_desired and not in_loop:
            structure = True
            continue
        elif residue_secstruct == ss_desired and not in_loop and structure:
            if i >= 2:
                start_position = i - 2 # include structured residues likely part of the beta bulge
            elif i == 1:
                start_position = i - 1
            in_loop = True
            structure = False
        elif residue_secstruct == ss_desired and in_loop:
            end_position = i + 1 # include structured residues likely part of the beta bulge
        elif residue_secstruct != ss_desired and in_loop:
            if end_position - start_position >= min_length:
                if end_position - start_position <= max_length or max_length == -1: #this maybe should be a single if statement
                    loops.append((start_position,end_position))
            in_loop = False;
            structure = True
    return loops

#if you only want 
def filter_structure_array()

#make a single list of all the residues
#for single heatmap
def get_all_res(loops):
    all_res = {}
    for loopset in loops:
        for i in range(loopset[0], loopset[1] + 1):
            all_res[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
    return all_res

#iterate through fasta file and count residues in loops
#for single heatmap
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
        all_res[loopset[1] + 2] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}

def get_dssp(pdb_file):
    import pyrosetta 
    pyrosetta.init()
    pose = pyrosetta.pose_from_file(pdb_file)
    struct_from_pose = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
    dssp_secstruct =struct_from_pose.get_dssp_secstruct()
    print(dssp_secstruct)
    return dssp_secstruct

def get_strand_data(pdb_file):
    import pyrosetta
    from pyrosetta.rosetta.core.select import get_residues_from_subset
    pyrosetta.init()
    pose = pyrosetta.pose_from_file(pdb_file)
    struct_from_pose = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
    dssp_secstruct = struct_from_pose.get_dssp_secstruct()
    strands = get_structure(dssp_secstruct, 'E', min_length = 7)
    print len(strands)
    core_selector = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    core_selector.set_layers(True, True, False)
    core_list = get_residues_from_subset(core_selector.apply(pose))
    print(core_list)
    print(dssp_secstruct)

def make_empty_dictionaries():
    two_res_loops = {}
    three_res_loops = {}
    four_res_loops = {}
    for i in range(-2, 5):
        if i < 3:
            two_res_loops[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
            three_res_loops[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
            four_res_loops[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
        elif i < 4:
            three_res_loops[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
            four_res_loops[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
        else:
            four_res_loops[i] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
    dicts = (two_res_loops, three_res_loops, four_res_loops)
    return dicts 

def check_loop_length(start,end,loops):
    if end - start == 1: #two residue loops
        loops[0].append((start,end))
    if end - start == 2: #three residue loops
        loops[1].append((start,end))
    if end - start == 3: #four residue loops
        loops[2].append((start,end))

#search secondary structure 
def sort_loops(dssp_secstruct):
    in_loop = False
    structure = False
    start_position = 0
    end_position = 0
    loops = [[],[],[]]
    for i, residue_secstruct in enumerate(dssp_secstruct):
        if residue_secstruct != 'L' and not in_loop:
            structure = True
            continue
        elif residue_secstruct == 'L' and not in_loop and structure:
            start_position = i 
            in_loop = True
            structure = False
        elif residue_secstruct == 'L' and in_loop:
            end_position = i  # include structured residues likely part of the beta bulge
        elif residue_secstruct != 'L' and in_loop:
            check_loop_length(start_position,end_position, loops)
            in_loop = False
            structure = True
    return loops

def compare_loop_position(fasta_name, loops, dicts):
    fasta = open(fasta_name, "r")
    for line in fasta.readlines():
        if line != [] and line[0] != '>':
            for i,loop_length in enumerate(loops): #iterates through 2,3, and 4 residue loopsets
                for loopset in loop_length: #gets each (start,end) set in 
                    position = -2
                    if loopset[0] - 2 >= 0:
                        for j in range(loopset[0] - 2, loopset[1] + 2):
                            if line[j] in dicts[i][position]:
                                dicts[i][position][line[j]] += 1
                            position += 1


def make_plots(dicts):
    for i, residue_count_dictionary in enumerate(dicts):
        name = str(i+2) + "_residue_loop"
        make_dataframe(residue_count_dictionary, name)

def make_dataframe(reside_counts,fasta_name):
    # dataframe = pd.DataFrame(reside_counts)
    # heatmap = sns.heatmap(dataframe, cbar = False, robust = True, annot = True)
    # #heatmap = sns.heatmap(dataframe, cmap="YlGnBu")
    # #heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=90, xticklabels = 2)
    # fig = heatmap.get_figure()
    # plot_name = clean_file_name(fasta_name) + ".png"
    # fig.savefig(plot_name)

    plt.figure()
    dataframe = pd.DataFrame(reside_counts)
    heatmap = sns.heatmap(dataframe, robust = True)
    plot_name = clean_file_name(fasta_name) + ".png"
    plt.savefig(plot_name)
     
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Takes in a list of pdb and fasta alignment files, uses the pdb secondary structure to get 2,3, and 4 resiue loop information from the alignment files.")
    parser.add_argument('FILES', type=str, nargs=1, help="file with sets of 'pdb_file_name fasta_file_name' on each line")
    parser.add_argument('-p', type=str, nargs=1, help="pdbs to be searched through. Can be /path/to/files/*.pdb or /path/to/your_file.pdb or your.pdb etc")
    parser.add_argument('-f', type=str, nargs=1, help="fasta file of aligned proteins to compare")
    parser.add_argument('-o','--out', type=str, nargs=1, default=[""], help="prefix for the output file")
                    
    args = parser.parse_args()
    input_file_name = args.FILES[0]
    if args.p != None: pdb = args.p[0] #name of PDB to extract the sequence from
    tag = args.out[0] #name of the output pdb

    #make alignment files of 2,3, and 4 residue loops
    dicts = make_empty_dictionaries()
    loop_profile = False
    strand_profile = True
    strand_placement_profile = False

    if(loop_profile):
        input_file = open(input_file_name, "r")
        for line in input_file.readlines():
            pdb_fasta = re.split(" |\n", line) #split line into pdb_name and fasta_name
            print(pdb_fasta)
            dssp_string = get_dssp(pdb_fasta[0]) #get string of secondary structure from pdb
            #dssp_string = pdb_fasta[0]
            loops = sort_loops(dssp_string)     #use 2nd struct string to determine whhere 2, 3, and 4 res loops
            print(loops)
            compare_loop_position(pdb_fasta[1], loops, dicts)
            make_plots(dicts)

    if (strand_profile):
        get_strand_data(input_file_name)
    # for i, pdb_file in enumerate(glob.glob(pdb)):
    #     dssp_string = get_dssp(pdb_file)

    #     #single alignment file
    #     loops = get_structure(dssp_string, "L")
    #     print (loops)
    #     all_res = get_all_res(loops)
    #     compare(args.Fasta[0],all_res, loops)
    #     make_dataframe(all_res, args.f[0])
