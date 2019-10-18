import re
import sys
import numpy as np
import networkx as nx
import os

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

def PBC3DF(c1, c2):
    """
        c1 and c2 are coordinates, either numpy arrays or lists
    """
    diffa = c1[0] - c2[0]
    diffb = c1[1] - c2[1]
    diffc = c1[2] - c2[2]

    if diffa > 0.5:
        c2[0] = c2[0] + 1.0
    elif diffa < -0.5:
        c2[0] = c2[0] - 1.0
    
    if diffb > 0.5:
        c2[1] = c2[1] + 1.0
    elif diffb < -0.5:
        c2[1] = c2[1] - 1.0
 
    if diffc > 0.5:
        c2[2] = c2[2] + 1.0
    elif diffc < -0.5:
        c2[2] = c2[2] - 1.0
    
    return c2

def bbelems(cifname, direc):

	path = os.path.join(direc, cifname)

	with open(path, 'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	elems = []
	elems_append = elems.append
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
			elems_append(s[1])

	return elems

def bb2array(cifname, direc):

	path = os.path.join(direc, cifname)

	with open(path, 'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	fcoords = []
	fcoords_append = fcoords.append
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
			fvec = np.array(map(float, s[2:5]))
			fcoords_append([s[0],fvec])

	pi = np.pi
	a,b,c,alpha,beta,gamma = map(float, (a,b,c,alpha,beta,gamma))
	ax = a
	ay = 0.0
	az = 0.0
	bx = b * np.cos(gamma * pi / 180.0)
	by = b * np.sin(gamma * pi / 180.0)
	bz = 0.0
	cx = c * np.cos(beta * pi / 180.0)
	cy = (c * b * np.cos(alpha * pi /180.0) - bx * cx) / by
	cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
	unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T
	
	norm_atom, norm_vec = fcoords[0]
	
	ccoords = [[n[0],np.dot(unit_cell, PBC3DF(norm_vec, n[1]))] for n in fcoords]
	com = np.average(np.array([n[1] for n in ccoords if re.sub('[0-9]','',n[0]) == 'X']), axis = 0)
	sccoords = [[n[0], n[1] - com] for n in ccoords]

	return sccoords

def bbbonds(cifname, direc):

	path = os.path.join(direc, cifname)

	with open(path, 'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	bonds = []
	bonds_append = bonds.append
	for line in cif:
		s = line.split()
		if isbond(s):
			bonds_append(s)
			
	return bonds

def X_vecs(cifname, direc, label):

	path = os.path.join(direc, cifname)

	with open(path, 'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	fcoords = []
	fcoords_append = fcoords.append
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
		if iscoord(s) and 'X' in s[0]:
			fcoords_append([s[0],np.array(map(float, s[2:5]))])

	pi = np.pi
	a,b,c,alpha,beta,gamma = map(float, (a,b,c,alpha,beta,gamma))
	ax = a
	ay = 0.0
	az = 0.0
	bx = b * np.cos(gamma * pi / 180.0)
	by = b * np.sin(gamma * pi / 180.0)
	bz = 0.0
	cx = c * np.cos(beta * pi / 180.0)
	cy = (c * b * np.cos(alpha * pi /180.0) - bx * cx) / by
	cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
	unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T

	mic_fcoords = [[vec[0],PBC3DF(fcoords[0][1],vec[1])] for vec in fcoords]

	if label:
		ccoords = [[vec[0],np.dot(unit_cell,vec[1])] for vec in mic_fcoords]
		com = np.average(np.asarray([vec[1] for vec in ccoords]), axis=0)
		shifted_ccoords = [[vec[0],vec[1] - com] for vec in ccoords]
	else:
		ccoords = [np.dot(unit_cell,vec[1]) for vec in mic_fcoords]
		com = np.average(ccoords, axis=0)
		shifted_ccoords = [vec - com for vec in ccoords]

	return shifted_ccoords

def bbcharges(cifname, direc):

	path = os.path.join(direc, cifname)

	with open(path, 'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	charges = []
	charges_append = charges.append
	elements = []
	elements_append = elements.append
	for line in cif:
		s = line.split()
		if iscoord(s):
			charges_append(s[-1])
			elements_append(s[1])
				
	return charges, elements

def calc_edge_len(cifname, direc):

	path = os.path.join(direc, cifname)

	with open(path, 'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	fcoords = []
	fcoords_append = fcoords.append
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
		if iscoord(s) and 'X' in s[0]:
			fcoords_append([s[0],np.array(map(float, s[2:5]))])

	pi = np.pi
	a,b,c,alpha,beta,gamma = map(float, (a,b,c,alpha,beta,gamma))
	ax = a
	ay = 0.0
	az = 0.0
	bx = b * np.cos(gamma * pi / 180.0)
	by = b * np.sin(gamma * pi / 180.0)
	bz = 0.0
	cx = c * np.cos(beta * pi / 180.0)
	cy = (c * b * np.cos(alpha * pi /180.0) - bx * cx) / by
	cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
	unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T

	mic_fcoords = [PBC3DF(fcoords[0][1],vec[1]) for vec in fcoords]
	ccoords = [np.dot(unit_cell,vec) for vec in mic_fcoords]

	if len(ccoords) > 2:
		print 'The edge cif', cifname, 'has to many connection points (Xs)'
		print 'Exiting'
		sys.exit()

	return np.linalg.norm(ccoords[0] - ccoords[1])

def cncalc(cifname, direc, cn1):

	path = os.path.join(direc, cifname)

	with open(path, 'r') as cif:
		cif = cif.read()
		cif = filter(None, cif.split('\n'))
	
	cn = 0
	nc = 0
	for line in cif:
		s = line.split()
		if iscoord(s):
			nc += 1
			if re.sub('[^a-zA-Z]','',s[0]) == 'X':
				cn += 1
	return cn

	