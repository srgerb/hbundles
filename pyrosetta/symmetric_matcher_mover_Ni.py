from collections import namedtuple
import math
import statistics
import os
import argparse

import pyrosetta

from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.numeric import xyzVector_double_t, angle_degrees
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.protein_interface_design.movers import TryRotamers
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.chemical import ChemicalManager
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.kinematics import Stub, Jump
from pyrosetta.rosetta.numeric import rotation_matrix

class CSTDefinition(namedtuple('CSTDefinition', ['residue', 'ligand', 'atoms', 'constraints', 'check_atoms'])):
    """Define jump defition for residue ligand interaction/sampling via matcher like constraints.

    residue - 3 letter code for residue type to sample
    ligand - 3 letter code for ligand type to sample
    atoms - List of 6 atoms for defining jump from constraints [A1,A2,A3,B1,B2,B3]... Jump is from A3 to B1
    constraints - Lists of 6 constraints for defining jumps/sampling.
    check_atoms - List of 2 atoms in ligand or residue if no ligand to check with axis.

    Constraints example:
    [[dist(A3,B1), angle_degrees(A2,A3,B1), angle_degrees(A3,B1,B2),
      dihedral_degrees(A1,A2,A3,B1), dihedral_degrees(A2,A3,B1,B2), dihedral_degrees(A3,B1,B2,B3)],
     [dist(A3,B1), angle_degrees(A2,A3,B1), angle_degrees(A3,B1,B2),
      dihedral_degrees(A1,A2,A3,B1), dihedral_degrees(A2,A3,B1,B2), dihedral_degrees(A3,B1,B2,B3)]...]

    Full example:
    CSTDefinition(residue='HID',
                  ligand='ZNO',
                  atoms=['ND1', 'CE1', 'NE2', 'ZN', 'O1', 'H1'],
                  constraints=[[2.2, 120, 109.5, 180, dihedral, 180] for dihedral in range(0, 360, 10)],
                  check_atoms=['ZN', 'O1'])
    """

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        if input_file[i] == ".":
            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    return tag

def get_rotamers_for_res_in_pose(in_pose, target_residue, sfx, ex1=True, ex2=True, ex3=False, ex4=False, check_clashes=True):
        packTask=pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(in_pose)
        packTask.set_bump_check(check_clashes)
        resTask=packTask.nonconst_residue_task(target_residue)
        resTask.or_ex1( ex1 )
        resTask.or_ex2( ex2 )
        resTask.or_ex3( ex3 )
        resTask.or_ex4( ex4 )
        resTask.restrict_to_repacking();
        packer_graph=pyrosetta.rosetta.utility.graph.Graph( in_pose.size() ) 
        rsf=pyrosetta.rosetta.core.pack.rotamer_set.RotamerSetFactory()
        rotset_ = rsf.create_rotamer_set( in_pose)
        rotset_.set_resid(target_residue)
        rotset_.build_rotamers(in_pose,
                           sfx,
                           packTask,
                           packer_graph,
                           False )
        print ("Found num_rotamers:", rotset_.num_rotamers())
        rotamer_set=[]
        for irot in range(1,rotset_.num_rotamers()+1):
            rotamer_set.append(rotset_.rotamer(irot).clone())
        return rotamer_set

class SymmetricMatcherMover():
    def __init__(self, cst_definitions, residue_selector=TrueResidueSelector(),
                 check_axis=[0, 0, 1], thresholds=[10.0, 0.5], symmetry_mover=None,
                 output_path='output/'):
        self.cst_definitions = cst_definitions
        self.residue_selector = residue_selector
        self.check_axis = xyzVector_double_t(*check_axis)
        self.thresholds = thresholds
        self.symmetry_mover = symmetry_mover
        self.output_path = output_path

    def _scan_rotamers(self, pose):
        rotamer_set = get_rotamers_for_res_in_pose(in_pose=pose, target_residue=self.residue_index,
                                                   sfx=pyrosetta.get_score_function(), check_clashes=False)
        for rotamer in rotamer_set:
            pose_clone = pose.clone()
            pose_clone.replace_residue(self.residue_index, rotamer, True)
            yield pose_clone

    def _scan_ligand(self, pose):
        cm = ChemicalManager.get_instance()
        rts = cm.residue_type_set(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)
        ligand = ResidueFactory.create_residue(rts.name_map(self.cst_definition.ligand))

        pose.append_residue_by_jump(ligand, self.residue_index,
                                     self.cst_definition.atoms[2],
                                     self.cst_definition.atoms[3])

        stubA = Stub(pose[self.residue_index].atom(self.cst_definition.atoms[2]).xyz(),
                     pose[self.residue_index].atom(self.cst_definition.atoms[1]).xyz(),
                     pose[self.residue_index].atom(self.cst_definition.atoms[0]).xyz())

        for constraints in self.cst_definition.constraints:
            # Logic from Core::Jump::from_bond_cst
            b1 = stubA.spherical(math.radians(constraints[3]), # /*dihedralA*/
                                 math.radians(180 - constraints[1]), # /*angleA*/
                                 constraints[0]) # /*disAB*/

            stubB1 = Stub(b1, pose[self.residue_index].atom(self.cst_definition.atoms[2]).xyz(), pose[self.residue_index].atom(self.cst_definition.atoms[1]).xyz())

            b2 = stubB1.spherical(math.radians(constraints[4]), # /*dihedralAB*/
                                  math.radians(180 - constraints[2]), # /*angleB*/
                                  pose[-1].atom(self.cst_definition.atoms[3]).xyz().distance(pose[-1].atom(self.cst_definition.atoms[4]).xyz()))

            stubB2 = Stub(b2, b1, pose[self.residue_index].atom(self.cst_definition.atoms[2]).xyz())

            b3 = stubB2.spherical(math.radians(constraints[5]), # /*dihedralB*/
                                  pyrosetta.rosetta.numeric.angle_radians(pose[-1].atom(self.cst_definition.atoms[3]).xyz(), pose[-1].atom(self.cst_definition.atoms[4]).xyz(), pose[-1].atom(self.cst_definition.atoms[5]).xyz()),
                                  pose[-1].atom(self.cst_definition.atoms[4]).xyz().distance(pose[-1].atom(self.cst_definition.atoms[5]).xyz()))

            stubB = Stub(b1, b2, b3)
            self.jump = Jump(stubA, stubB)
            pose.set_jump(1, self.jump)
            yield pose

    def _check(self, pose, ligand=True):
        if ligand:
            orientation = pose[-1].atom(self.cst_definition.check_atoms[0]).xyz().z - pose[-1].atom(self.cst_definition.check_atoms[1]).xyz().z
            if orientation < 0:
                 return False
            angle_error = min(angle_degrees(xyzVector_double_t(0,0,0), self.check_axis,
                                            pose[-1].atom(self.cst_definition.check_atoms[0]).xyz(),
                                            pose[-1].atom(self.cst_definition.check_atoms[1]).xyz()),
                              abs(180 - angle_degrees(xyzVector_double_t(0,0,0), self.check_axis,
                                                      pose[-1].atom(self.cst_definition.check_atoms[0]).xyz(),
                                                      pose[-1].atom(self.cst_definition.check_atoms[1]).xyz())))

            distance_error = statistics.mean([xyzVector_double_t(0,0,pose[-1].atom(self.cst_definition.check_atoms[0]).xyz().z).distance_squared(pose[-1].atom(self.cst_definition.check_atoms[0]).xyz()),
                                             xyzVector_double_t(0,0,pose[-1].atom(self.cst_definition.check_atoms[1]).xyz().z).distance_squared(pose[-1].atom(self.cst_definition.check_atoms[1]).xyz())])

        else:
            angle_error = min(angle_degrees(xyzVector_double_t(0,0,0), self.check_axis,
                                            pose[self.residue_index].atom(self.cst_definition.check_atoms[0]).xyz(),
                                            pose[self.residue_index].atom(self.cst_definition.check_atoms[1]).xyz()),
                              abs(180 - angle_degrees(xyzVector_double_t(0,0,0), self.check_axis,
                                                      pose[self.residue_index].atom(self.cst_definition.check_atoms[0]).xyz(),
                                                      pose[self.residue_index].atom(self.cst_definition.check_atoms[1]).xyz())))

            distance_error = statistics.mean([xyzVector_double_t(0,0,pose[self.residue_index].atom(self.cst_definition.check_atoms[0]).xyz().z).distance_squared(pose[self.residue_index].atom(self.cst_definition.check_atoms[0]).xyz()),
                                             xyzVector_double_t(0,0,pose[self.residue_index].atom(self.cst_definition.check_atoms[1]).xyz().z).distance_squared(pose[self.residue_index].atom(self.cst_definition.check_atoms[1]).xyz())])

        errors = [angle_error, distance_error]
        return all([error < threshold for error, threshold in zip(errors, self.thresholds)])

    def apply(self, pose, pdb_name):
        for i, cst_definition in enumerate(self.cst_definitions):
            self.cst_definition = cst_definition
            residue_list = get_residues_from_subset(self.residue_selector.apply(pose))

            for j in get_residues_from_subset(self.residue_selector.apply(pose)):
                self.residue_index = j
                dumped = False

                pose_clone1 = pose.clone()
                MutateResidue(self.residue_index, self.cst_definition.residue).apply(pose_clone1)

                for k, pose_clone2 in enumerate(self._scan_rotamers(pose_clone1)):
                    if dumped == True:
                        break
                    if self.cst_definition.ligand:
                        for l, pose_clone3 in enumerate(self._scan_ligand(pose_clone2)):
                            if self._check(pose_clone3):
                                pose_clone3.dump_pdb(os.path.join(self.output_path, '_'.join([pdb_name, str(i), str(j), str(k), str(l), 'asymm']) + '.pdb'))
                                if self.symmetry_mover:
                                    pose_clone4 = pose_clone3.clone()
                                    self.symmetry_mover.apply(pose_clone4)
                                    pose_clone4.dump_pdb(os.path.join(self.output_path, '_'.join([pdb_name, str(i), str(j), str(k), str(l), 'symm']) + '.pdb'))
                                    dumped = True
                                    break

                    else:
                        if self._check(pose_clone2, ligand=False):
                            pose_clone2.dump_pdb(os.path.join(self.output_path, '_'.join([pdb_name, str(i), str(j), str(k), 'asymm']) + '.pdb'))
                            if self.symmetry_mover:
                                pose_clone3 = pose_clone2.clone()
                                self.symmetry_mover.apply(pose_clone3)
                                pose_clone3.dump_pdb(os.path.join(self.output_path, '_'.join([pdb_name, str(i), str(j), str(k), 'symm']) + '.pdb'))

if __name__ == '__main__':
    import glob

    parser = argparse.ArgumentParser(description="Takes in any number of pdbs and searches through to find a symmetric metal site based on constraints in this file.")
    parser.add_argument('PDBS', type=str, nargs=1, help="pdbs to be searched through. Can be /path/to/files/*.pdb or /path/to/your_file.pdb or your.pdb etc")
    parser.add_argument('-o','--out', type=str, nargs=1, default=[""], help="prefix for the output file")

    args = parser.parse_args()
    input_files=args.PDBS[0] #name of PDB to extract the sequence from
    tag = args.out[0] #name of the output pdb
    pyrosetta.init('-extra_res_fa /home/srgerb/hBundles/nickel/Ni_C4V.params /home/srgerb/hBundles/nickel/HID.params -out:file:output_virtual true')

    # More CSTDefinitions could be added to the list...
    #cst_definitions = [CSTDefinition(residue='HID',
    #                                 ligand='NIV',
    #                                 atoms=['CG', 'CD2', 'NE2', 'Ni', 'V1', 'V2'],
    #                                 constraints=[[2.2, 120, 105, 180.0, dihedral, 180] for dihedral in range(0, 360, 10)],
    #                                 check_atoms=['Ni', 'V1']),]    
    # For example:
    cst_definitions =[  CSTDefinition(residue='HID',
                                     ligand='NIV',
                                     atoms=['CG', 'CD2', 'NE2', 'Ni', 'V1', 'V2'],
                                     constraints=[[2.2, 120, 105, 180.0, dihedral, 180] for dihedral in range(0, 360, 10)],
                                     check_atoms=['Ni', 'V1'])] # Check the same atoms
    #surface_selector = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    #surface_selector.set_layers(False, False, True)
 
    symmetry_mover = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover('/home/srgerb/hBundles/SymFiles/C4_Z.sym')
    smm = SymmetricMatcherMover(cst_definitions=cst_definitions,
                                symmetry_mover=symmetry_mover,
                                output_path=tag,
                                thresholds=[15, 0.5])

    for i, pdb_file in enumerate(glob.glob(input_files)):
        pose = pyrosetta.pose_from_file(pdb_file)
        smm.apply(pose, clean_file_name(pdb_file))
