#usr/bin/python
#calculate rmsd between 2 sets of coordinates.

from math import sqrt
import sys
import os

parsing = False
coords = {}
d = {}
dd = {}
tbl = {}
conv_GvA_ls = []
euc1 = 0
euc_tot = 0
rmsd = 0
no_atoms = 0
count = 1
out_dir = "/Users/nj001/Documents/Docking/201604/20160425/flex_analysis/open_3/"
ATOM_conv_dir = ("/Users/nj001/Documents/Docking/ATOM_conversion_files/methyl_substituted_compounds_noH/")
symm = ['4MP','26MP','35MP','2356MP','23456MP']

os.chdir(ATOM_conv_dir)
conversion = [f for f in os.listdir(".") if f.startswith("ATOM_conversion_")]
os.chdir(out_dir)
GvAs = [f for f in os.listdir(".") if f.endswith("GvA.mol2")]
conversion = sorted(conversion)
GvAs = sorted(GvAs)
conv_GvA_ls.append(GvAs)
conv_GvA_ls.append(conversion)
print conv_GvA_ls
symmetrical_ph = False

for x in range(len(GvAs)):
	
	GvA = conv_GvA_ls[0][x]
	methyl_subs = GvA[:-22]
	
	if methyl_subs in symm:
		symmetrical_ph = True
	else:
		symmetrical_ph = False
	
	conv = conv_GvA_ls[1][x]
	f = open(GvA,'r')  #coordinate file to parse (mol2 format)
	atom_conv = open(ATOM_conv_dir+conv,'r')  #file listing which atoms are which for vina vs GOLD output
	
	#### PARSE FILE FOR ATOM COORDINATES ####
	for line in f:
	
		#read only the relevant lines in the file (parsing between the 
		#start of atom list and start of bond list)
		if line.startswith("@<TRIPOS>BOND"):
			parsing = False
			#make a dictionary of the dictionaries of coordinates
			dd[count] = d
			d = {}
			count += 1
		if parsing:
			ln = line.split()
			#add the coordinates to a dictonary with the atom number as the key
			atom = int(ln[0])
			coords['atom_type'] = ln[1][0]
			coords['x']=float(ln[2])
			coords['y']=float(ln[3])
			coords['z']=float(ln[4])
			d[atom]=coords
			coords={}
		if line.startswith("@<TRIPOS>ATOM"):
			parsing = True

	f.close()
	no_vina_sols = count - 50
	count = 1

	#### USE ATOM CONVERSION TO CALCULATE RMSD ####

	for line in atom_conv: #create dictionary to look up which atoms relate to which 
		if not line.startswith('#'):
			no_atoms += 1
			v=line.split()
			tbl[int(v[1])]=int(v[0])
	
	print 'Gva = ',GvA, 'conv = ',conv,'no_atoms = ', no_atoms	
		
	with open('rmsds_'+GvA[:-9]+'.txt','w') as outfile:
		
		header1 = ("#",GvA[:-9],'\n')
		header2 = ('#VINA GOLD_chimera GOLD_soln RMSD(A)\n')
		
		outfile.writelines(header1)
		outfile.writelines(header2)
		
		for i in range(1,(no_vina_sols)): #iterate vina results
			for j in range((no_vina_sols),(no_vina_sols+50)): #iterate gold results
				for key in tbl: #iterate vina result atoms
					b = tbl[key]
					a = key
					
					if not dd[i][a]['atom_type'] == 'H':
						for l in ['x','y','z']:
							euc1 += (dd[i][a][l]-dd[j][b][l])*(dd[i][a][l]-dd[j][b][l])
						euc_tot += euc1
						euc1 = 0
				rmsd = sqrt(euc_tot/(no_atoms))
				line = (str(i),' ',str(j),' ',str(j-no_vina_sols+1),' ',str(rmsd),'\n')
				outfile.writelines(line)
				rmsd = 0
				euc_tot = 0
	no_atoms = 0
