import pyrosetta
import argparse
from pyrosetta.rosetta.core.pose import make_pose_from_sequence

parser = argparse.ArgumentParser(description="takes in a .params file and spits out a pdb")
parser.add_argument('PARAM', type=str, nargs=1, help="param file to turn into pdb")
parser.add_argument('-o', type=str, nargs=1, default= ["test.pdb"], help= "name of output file")
args = parser.parse_args()

init= '-extra_res_fa ' + args.PARAM[0] + ' -out:file:output_virtual true'
pyrosetta.init(init)
pose = pyrosetta.Pose()
make_pose_from_sequence(pose, sequence='X[FE3]', type_set_name="fa_standard", auto_termini=False)
pose.dump_pdb(args.o[0])
