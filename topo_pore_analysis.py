import glob
import re
import numpy as np
import itertools

from scipy.spatial import Voronoi
from ase.io import read

topos = ['tty' , 'hwx', 'ucp' , 'sab', 'sqca', 
		 'edq' , 'apo', 'lnj' , 'scu', 'qza' , 
		 'pts' , 'hof', 'cut' , 'ftw', 'sit' , 
		 'lil' , 'lim', 'svnd', 'phx', 'pcu' , 
		 'ddi' , 'fog', 'sur' , 'rtl', 'csq' , 
		 'rhrb', 'reo', 'lvtb', 'mcn', 'dag' , 
		 'fcu' , 'tfk', 'pyr' , 'pfm', 'rht' , 
		 'tfb' , 'kkm', 'nbob', 'ceq', 'flu' , 
		 'ssa' , 'fsc', 'stx' , 'xbq', 'ibe' , 
		 'acs' , 'stu', 'xll' , 'stp', 'xly' , 
		 'act' , 'bcu', 'the' , 'xbf', 'tru' , 
		 'sty' , 'mrc']

templates = ['svnd.cif']
templates = [c for c in glob.glob('*.cif') if c.split('.')[0] in topos]

def PBC3DF_sym(vec1, vec2):

	dX,dY,dZ = vec1 - vec2
			
	if dX > 0.5:
		s1 = 1
		ndX = dX - 1.0
	elif dX < -0.5:
		s1 = -1
		ndX = dX + 1.0
	else:
		s1 = 0
		ndX = dX
				
	if dY > 0.5:
		s2 = 1
		ndY = dY - 1.0
	elif dY < -0.5:
		s2 = -1
		ndY = dY + 1.0
	else:
		s2 = 0
		ndY = dY
	
	if dZ > 0.5:
		s3 = 1
		ndZ = dZ - 1.0
	elif dZ < -0.5:
		s3 = -1
		ndZ = dZ + 1.0
	else:
		s3 = 0
		ndZ = dZ

	return np.array([ndX,ndY,ndZ]), np.array([s1,s2,s3])

def Voronoi_tessalate(atoms):

	base_coords = [a.position for a in atoms]
	repeat_unit_cell = atoms.get_cell().T
	inv_ruc = np.linalg.inv(repeat_unit_cell)
	basis = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]
	mesh = []

	for coord in base_coords:

		mesh.append(coord)
		
		fcoord = np.dot(inv_ruc, coord)
		zero_threshold_indices = fcoord < 1e-6
		fcoord[zero_threshold_indices] = 0.0
		one_threshold_indices = abs(fcoord - 1.0) < 1e-6
		fcoord[one_threshold_indices] = 0.0

		if np.all(fcoord):
			trans_vecs = [-1 * b for b in basis] + basis
		else:
			trans_vecs = [basis[dim] for dim in (0,1,2) if fcoord[dim] == 0.0]
			combs = list(itertools.combinations(trans_vecs, 2)) + list(itertools.combinations(trans_vecs, 3))
			for comb in combs:
				compound = np.array([0.0,0.0,0.0])
				for vec in comb:
					compound += vec
				trans_vecs.append(compound)

		for vec in trans_vecs:
			trans_coord = [np.round(i, 6) for i in np.dot(repeat_unit_cell, fcoord + vec)]
			mesh.append(trans_coord)

	mesh = np.asarray(mesh)
	vor = Voronoi(mesh)

	return vor, mesh

for cif in templates:
	
	template = read(cif)
	unit_cell = template.get_cell().T
	inv_uc = np.linalg.inv(unit_cell)
	
	coords = []
	labels = []

	for v in template:
		coords.append(v.position)
		labels.append(v.symbol)

	nv = len(set(labels))

	vor, mesh = Voronoi_tessalate(template)
	voronoi_vertices = vor.vertices

	in_cell_vertices = []
	count = 0
	for vertex in voronoi_vertices:
		fcoord = np.dot(inv_uc, vertex)
		
		include = True
		for dim in fcoord:
			if dim > 1 or dim < 0:
				include = False
		
		if include:
			in_cell_vertices.append(vertex)

	spectrums = []
	norm_val = 0
	for v in in_cell_vertices:
		
		dists = []
		vf = np.dot(inv_uc, v)

		for nv in coords:

			nvf = np.dot(inv_uc, nv)
			fdist, sym = PBC3DF_sym(vf, nvf)
			dist = np.dot(unit_cell, fdist)
			dists.append(np.linalg.norm(dist))
		
		md = min(dists)
		spectrum = [d for d in dists if (d - md) < 1e-5]
		spectrums.append(spectrum)

		if spectrum[0] > norm_val:
			norm_val = spectrum[0]


	spec_vals = [s[0] for s in spectrums]/norm_val

	avg_diff = 0
	for i in range(len(spec_vals)):
		ival = spec_vals[i]
		for j in range(i + 1, len(spec_vals)):
			jval = spec_vals[j]
			diff = abs(ival - jval)
			avg_diff += diff

	avg_diff /= float(len(spec_vals))
	num_nodes = [len(s) for s in spectrums]
	max_c = max(num_nodes)
	min_c = min(num_nodes)
	nsites = len(set(num_nodes))

	print cif.split('.')[0], np.round(avg_diff,5), nsites


	





