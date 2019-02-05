import numpy as np
import re

metals = ['Mg', 'Al', 'Ca', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo',
		  'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Ba', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 
		  'Pb', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U' ]

def fix_charges(placed_all):
	placed_all = np.asarray(placed_all)
	net_charge = np.sum(map(float, placed_all[:,4]))
	
	dont_change_charge_atoms = []
	dont_change_charge_atoms_append = dont_change_charge_atoms.append
	for l in placed_all:
		elem = re.sub('[0-9]','',l[0])
		if elem in metals:
			dont_change_charge_atoms_append(l[0])
		elif 'FG' in elem or 'fg' in elem or 'fG' in elem or 'Fg' in elem:
			dont_change_charge_atoms_append(l[0])

	num_atoms = float(len(placed_all) - len(dont_change_charge_atoms))
	rcb = net_charge/num_atoms

	fc_placed_all = []
	fc_placed_all_append = fc_placed_all.append
	for l in placed_all:
		e,x,y,z,c,oe,i = l
		c = float(c)
		if e in dont_change_charge_atoms:
			fc_placed_all_append([e,x,y,z,c,oe,i])
		else:
			nc = c - rcb
			fc_placed_all_append([e,x,y,z,nc,oe,i])

	nnet_charge = 0
	for l in fc_placed_all:
		nnet_charge = nnet_charge + l[4]

	return fc_placed_all, nnet_charge, net_charge, rcb
