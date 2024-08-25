#!/usr/bin/env python

import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np
import argparse
from tqdm import tqdm
from scipy.spatial.transform import Rotation


parser = argparse.ArgumentParser()
parser.add_argument('--pdb',     required=True,  type=str,
                                        help='path to pdb file')
parser.add_argument('--dcd',     required=True,  type=str,
                                        help='path to dcd file')

args = parser.parse_args()

# Load the universe with your topology and trajectory files
u = mda.Universe(args.pdb, args.dcd)

# Select the first and last carbon atoms in the polyethylene chain
first_atom = u.select_atoms('name C')[0]  # Assumes 'name C' selects all carbon atoms, change as necessary
last_atom = u.select_atoms('name C')[-1]  # Selects the last carbon atom


#"""
#Aligns the chain so that the first atom is at the origin and the vector
#from the first atom to the last atom lies along the x-axis.
#
#Parameters:
#- u: MDAnalysis Universe object
#- first_atom: Atom object representing the first carbon atom
#- last_atom: Atom object representing the last carbon atom
#"""

origin = first_atom.position
# Step 1: Translate the system so the first carbon atom is at the origin
trans = transformations.center_in_box(u.select_atoms(f'id {first_atom.id}'),
                                      point=[0, 0, 0])
u.trajectory.add_transformations(trans)

# Save the transformed trajectory to a new file
with mda.Writer('aligned.dcd', u.atoms.n_atoms) as W:
    for t in tqdm(u.trajectory):
        # Step 2: Calculate the vector from the first atom to the last atom
        vector = last_atom.position - first_atom.position
        # Normalize the vector
        direction = vector / np.linalg.norm(vector)
        y_axis = np.array([0, 1, 0])
    
        theta = np.arccos(np.dot(y_axis, direction))
        axis=np.cross(direction, y_axis)
        axis/=np.linalg.norm(axis)
        # Step 4: Rotate such that C1-Cend align with 0,1,0
        rot=Rotation.from_rotvec(theta*axis)
        u.atoms.positions = rot.apply(u.atoms.positions)
        W.write(u)
