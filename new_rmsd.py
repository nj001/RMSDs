# usr/bin/python
# calculate rmsd between 2 sets of coordinates which are for the same molecule
# but atom indexes are assigned differently for the same atoms
# to work with .pdb files with multiple structures in them
#### master_atoms[model][atom][x, y, z, atm_type] ###
#### master_connects[model][connection_list][atm_index, bonded_atm_1, ... bonded_atm_n] ###
# result [atms/connects][model][atm_idx][coords/bonded_atm_indexes]

import numpy as np
import re
import os
import sys
from math import sqrt

def get_HETATM(line, dictionary):
	'''
	take a line from a pdb file that is a HETATM record
	get coordinates and append to dictionary with atom 
	index as the key
	'''
	lin = line.split()
	atm_idx = int(lin[1])
	x = float(lin[4])
	y = float(lin[5])
	z = float(lin[6])
	atm_type = lin[2]
	dictionary[atm_idx] = [x,y,z,atm_type])
	
def get_CONECT(line,dictionary):
	'''
	take a line from a pdb file that is a CONECT record
	get atom index and atoms bonded to it and append to
	list
	'''
	lin = line.split()
	atm_idx = int(lin[1])
	atms_bonded = []
	for i in range(2,len(lin)):
		atms_bonded.append(int(lin[i]))
	dictionary[atm_idx] = atms_bonded

def make_matrices(file):
	'''
	take a file and make list of dictionaries for atom 
	coords and bonds made by each atom
	return the 2 lists in a dictionary
	result = {'atoms', 'connects'}
	'''
	mol_atoms = {}
	mol_connects = {}
	master_atoms = []	
	master_connects = []
	for line in file:
		ln = line.split()
		if ln[0] == "HETATM":
			get_HETATM(line, mol_atoms)
		elif ln[0] == "ENDMDL":
			master_atoms.append(mol_atoms)
			mol_atoms = []
		elif ln[0] == "CONECT":
			get_CONECT(line, mol_connects)
		elif ln[0] == "MODEL" or ln[0] == "END":
			master_connects.append(mol_connects)
			mol_connects = []
	return {"atoms":master_atoms, "connects":master_connects}
	
def check_distance(model,atom1,atom2):
	tot = 0
	for i in range(3):
		tot += (model[atom1][i] - model[atom2][i]) ** 2
	euc_dist = sqrt(tot)
	return euc_dist
