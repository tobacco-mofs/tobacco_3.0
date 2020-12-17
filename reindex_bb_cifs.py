import re
import numpy as np
import os
import glob
from datetime import date

PT = ['H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
	  'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
	  'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe',
	  'Cs', 'Ba', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 
	  'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 
	  'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'FG', 'X' ]

def nn(string):
	return re.sub('[^a-zA-Z]','', string)

def nl(string):
	return re.sub('[^0-9]','', string)

def isfloat(value):
	"""
		determines if a value is a float
	"""
	try:
		float(value)
		return True
	except ValueError:
		return False

def iscoord(line):
	"""
		identifies coordinates in CIFs
	"""
	if nn(line[0]) in PT and line[1] in PT and False not in map(isfloat,line[2:5]):
		return True
	else:
		return False
	
def isbond(line):
	"""
		identifies bonding in cifs
	"""
	if nn(line[0]) in PT and nn(line[1]) in PT and isfloat(line[2]) and line[-1] in ('S', 'D', 'T', 'A'):
		return True
	else:
		return False

def reindex(direc, consider_charges=False):

	cifs = glob.glob(direc + os.sep + '*.cif')

	for path in cifs:

		coords_and_elems = []
		atom_names = []
		charges = []
		bonds = []
		filename = path.split(os.sep)[-1].split('.')[0]

		with open(path, 'r') as cif:

			cif = cif.read()
			cif = filter(None, cif.split('\n'))
		
		for line in cif:
			s = line.split()
			if '_cell_length_a' in line:
				a = s[1]
			if '_cell_length_b' in line:
				b = s[1]
			if '_cell_length_c' in line:
				c = s[1]
			if '_cell_angle_alpha' in line:
				alpha = s[1]
			if '_cell_angle_beta' in line:
				beta = s[1]
			if '_cell_angle_gamma' in line:
				gamma = s[1]
			
			if iscoord(s):
				
				coords_and_elems.append(s[1:])
				atom_names.append(s[0])

				if consider_charges:
					charges.append(s[-1])
				else:
					charges.append(0.0)

			if isbond(s):

				bonds.append(s)
		
		new_name_key = {}
		for i,name in enumerate(atom_names):
			new_name_key[name] = re.sub('[0-9]','',name) + str(i+1)
	
		with open(direc + os.sep + filename + '_reindexed.cif', 'w') as out:
	
			out.write('data' + '_' + filename + '\n')
			out.write('_audit_creation_date              ' + str(date.today()) + '\n')
			out.write("_audit_creation_method            Ryther's reindex function\n")
			out.write("_symmetry_space_group_name_H-M    'P1'\n")
			out.write('_symmetry_Int_Tables_number       1\n')
			out.write('_symmetry_cell_setting            triclinic\n')
			out.write('loop_\n')
			out.write('_symmetry_equiv_pos_as_xyz\n')
			out.write('  x,y,z\n')
			out.write('_cell_length_a                    ' + a + '\n')
			out.write('_cell_length_b                    ' + b + '\n')
			out.write('_cell_length_c                    ' + c + '\n')
			out.write('_cell_angle_alpha                 ' + alpha + '\n')
			out.write('_cell_angle_beta                  ' + beta + '\n')
			out.write('_cell_angle_gamma                 ' + gamma + '\n')
			out.write('loop_\n')
			out.write('_atom_site_label\n')
			out.write('_atom_site_type_symbol\n')
			out.write('_atom_site_fract_x\n')
			out.write('_atom_site_fract_y\n')
			out.write('_atom_site_fract_z\n')
			out.write('_atom_site_charge\n')
	
			for old_name,coord,charge in zip(atom_names,coords_and_elems,charges):
				
				line = [new_name_key[old_name]] + coord + [str(charge)]
				out.write(' '.join(line))
				out.write('\n')

			out.write('loop_\n')
			out.write('_geom_bond_atom_site_label_1\n')
			out.write('_geom_bond_atom_site_label_2\n')
			out.write('_geom_bond_distance\n')
			out.write('_geom_bond_site_symmetry_2\n')
			out.write('_ccdc_geom_bond_type\n')
	
			for bond in bonds:
	
				#bond = bond.split()
				new_bond = [new_name_key[bond[0]], new_name_key[bond[1]]] + bond[2:]
				out.write(' '.join(new_bond))
				out.write('\n')





reindex('edges')
