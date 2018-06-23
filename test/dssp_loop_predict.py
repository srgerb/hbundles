#import pyrosetta        #dssp
import glob             #inputfiles (maybe)
import argparse         # parse arguments from command line
import pandas as pd
import seaborn as sns

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        if input_file[i] == ".":
            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
            tag = ""
    return tag

# def get_dssp(pdb_file):
#     pyrosetta.init('-out:file:output_virtual true')
#     pose = pyrosetta.pose_from_file(pdb_file)
#     struct_from_pose = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
#     dssp_secstruct =struct_from_pose.get_dssp_secstruct()
#     print(dssp_secstruct)

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
            start_position = i - 2 # include structured residues likely part of the beta bulge
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
    amino_acids = {"A:0","C:0","D:0","E:0","F:0","G:0","H:0","I:0","J:0","K:0","L:0","M:0","N:0","P:0","Q:0","R:0","S:0","T:0","V:0","W:0","Y:0"}
    for loopset in loops:
        for i in range(loopset[0], loopset[1] + 1):
            all_res[i] = amino_acids
    return all_res

def compare(fasta_name,all_res):
    fasta = open(fasta_name, "r")
    for line in fasta.readlines():
        if line != [] and line[0] != '>':
            for i,residue in enumerate(line):
                if i in all_res and residue in all_res[i]:
                    print ("increment " + str(i) + " "+ residue)
                    res = all_res[i][residue]
                    res += 1
                    all_res[i][residue] = res


def make_dataframe(all_res):
    data = []
    res_numbers = all_res.keys()
    columns = all_res[res_numbers[0]].keys()
    rescount.append(columns)
    for residue in res_numbers:
        rescount.append(list(all_res[residue].values))
    amino_acids = {"A:[]","C:[]","D:[]","E:[]","F:[]","G:[]","H:[]","I:[]","J:[]","K:[]","L:[]","M:[]","N:[]","P:[]","Q:[]","R:[]","S:[]","T:[]","V:[]","W:[]","Y:[]"}
    dataframe = pd.DataFrame(columns = amino_acids, data= data)

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
        #dssp_string = get_dssp(pdb_file)
        dssp_string = "LLLLLEEEEEEEELLLLLEEEEEEEEEEEEEELLLEEEEEEEEEEEELLEEEEEEEEEEEEEELLEEEEEEEEEEEELLLLEEEEEEEEEEEELLEEEEEEEEEELLLLEELLLLLLLEEEELLEEEEEEEEEELLEEEEEEEEEEEEELLLLELLLLLLLLELLLEEEEEEEEEEEELLEEEEEEEEELHHHLEEEEEEEEEEEELLLLHHHHHLHHHHHHHHLHHHHHHLLLLLLLLLLEEEELLLLLLLLLLLLL"
        loops = get_loops(dssp_string)
        print (loops)
        all_res = get_all_res(loops)
        print ("all_res[11]")
        print (all_res[11])
        print ("all_res[11][A]")
        print (all_res[11]['A'])
        #print(all_res.values())
        compare(args.Fasta[0],all_res)
        #print(all_res.values())
        # make_heatmap(all_res)
