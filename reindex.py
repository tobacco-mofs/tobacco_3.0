import sys
import os
import re
import datetime
import numpy as np
import glob

PT = ['H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
	  'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
	  'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe',
	  'Cs', 'Ba', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 
	  'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 
	  'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'FG', 'X' ]

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
	if (re.sub('[^a-zA-Z]','', line[0]) in PT and
	line[1] in PT and
	isfloat(line[2]) and
	isfloat(line[3]) and
	isfloat(line[4])):
		return True
	else:
		return False
	
def isbond(line):
	"""
		identifies bonding in cifs
	"""
	if (re.sub('[^a-zA-Z]','', line[0]) in PT and
	re.sub('[^a-zA-Z]','', line[1]) in PT and
	isfloat(line[2]) and
	not isfloat(line[3]) and
	not isfloat(line[4])):
		return True
	else:
		return False

def reindex(pathfile,charges):

	file = pathfile.split(os.sep)[-1]
	path_list = pathfile.split(os.sep)[0:-1]

	path = ''
	c = 0
	for p in path_list:
		c += 1
		if c == 1:
			path = path + p
		else:
			path = path + os.sep + p

	with open(pathfile,'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	coords = []
	coords_append = coords.append
	bonds = []
	bonds_append = bonds.append

	cc = 0

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
			cc += 1
			coords_append(s + [cc])
		if isbond(s):
			bonds_append(s)

	ind_dict = dict((k[0],k[-1]) for k in coords)

	for l in coords:
		l[0] = re.sub('[0-9]','',l[0]) + str(l[-1])
		l.pop(-1)

	#print ind_dict

	for l in bonds:
		l[0] = re.sub('[0-9]','',l[0]) + str(ind_dict[l[0]])
		l[1] = re.sub('[0-9]','',l[1]) + str(ind_dict[l[1]])

	new_path = os.path.join(path, file + '_ri')

	with open(new_path, 'w') as out:
		out.write('data_' + file[0:-4] + '\n')
		out.write('_audit_creation_date              ' + datetime.datetime.today().strftime('%Y-%m-%d') + '\n')
		out.write("_audit_creation_method            'tobacco_3.0'" + '\n')
		out.write("_symmetry_space_group_name_H-M    'P1'" + '\n')
		out.write('_symmetry_Int_Tables_number       1' + '\n')
		out.write('_symmetry_cell_setting            triclinic' + '\n')
		out.write('loop_' + '\n')
		out.write('_symmetry_equiv_pos_as_xyz' + '\n')
		out.write('  x,y,z' + '\n')
		out.write('_cell_length_a                    ' + str(a) + '\n')
		out.write('_cell_length_b                    ' + str(b) + '\n')
		out.write('_cell_length_c                    ' + str(c) + '\n')
		out.write('_cell_angle_alpha                 ' + str(alpha) + '\n')
		out.write('_cell_angle_beta                  ' + str(beta) + '\n')
		out.write('_cell_angle_gamma                 ' + str(gamma) + '\n')
		out.write('loop_' + '\n')
		out.write('_atom_site_label' + '\n')
		out.write('_atom_site_type_symbol' + '\n')
		out.write('_atom_site_fract_x' + '\n')
		out.write('_atom_site_fract_y' + '\n')
		out.write('_atom_site_fract_z' + '\n')
		out.write('_atom_site_U_iso_or_equiv' + '\n')
		out.write('_atom_site_adp_type' + '\n')
		out.write('_atom_site_occupancy' + '\n')
		
		if charges:
			out.write('_atom_site_charge' + '\n')

		for l in coords:
			name = l[0]
			elem = l[1]
			vec = map(float,l[2:5])

			if charges:
				charge = l[-1]
			else:
				charge = ''

			print charge

			out.write('{:7} {:>4} {:>15} {:>15} {:>15} {:>22} {:>15}'.format(name, elem, "%.10f" % np.round(vec[0],10), "%.10f" % np.round(vec[1],10), "%.10f" % np.round(vec[2],10),'0.00000  Uiso   1.00',charge))
			out.write('\n')

		out.write('loop_' + '\n')
		out.write('_geom_bond_atom_site_label_1' + '\n')
		out.write('_geom_bond_atom_site_label_2' + '\n')
		out.write('_geom_bond_distance' + '\n')
		out.write('_geom_bond_site_symmetry_2' + '\n')
		out.write('_ccdc_geom_bond_type' + '\n')

		for l in bonds:
			out.write('{:7} {:>7} {:>5} {:>7} {:>3}'.format(l[0], l[1], "%.3f" % float(l[2]), l[3], l[4]))
			out.write('\n')

def apply_reindex(CHARGES):

	ecifs = glob.glob('edges' + os.sep + '*.cif')
	ncifs = glob.glob('nodes' + os.sep + '*.cif')
	
	for e in ecifs:
		reindex(e,CHARGES)
	for n in ncifs:
		reindex(n,CHARGES)
	
	for e in ecifs:
		os.remove(e)
	for n in ncifs:
		os.remove(n)
	
	ri_ecifs = glob.glob('edges' + os.sep + '*.cif_ri')
	ri_ncifs = glob.glob('nodes' + os.sep + '*.cif_ri')
	
	for e in ri_ecifs:
		os.rename(e,e.split('.')[0]+'.cif')
	for n in ri_ncifs:
		os.rename(n,n.split('.')[0]+'.cif')

