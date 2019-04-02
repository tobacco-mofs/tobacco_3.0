import numpy as np
import networkx as nx
import re
import os
from ciftemplate2graph import isvert, PBC3DF

def omega2coords(augTG, sc_omega_plus, uc_params, num_vertices, template, g, CHECK):

	pi = np.pi
	sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = uc_params
	ax = sc_a
	ay = 0.0
	az = 0.0
	bx = sc_b * np.cos(sc_gamma * pi/180.0)
	by = sc_b * np.sin(sc_gamma * pi/180.0)
	bz = 0.0
	cx = sc_c * np.cos(sc_beta * pi/180.0)
	cy = (sc_c * sc_b * np.cos(sc_alpha * pi/180.0) - bx * cx) / by
	cz = (sc_c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
	sc_unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T

	path = os.path.join('templates', template)

	with open(path, 'r') as cif:

		cif = cif.read()
		cif = filter(None, cif.split('\n'))

	discard_coords = []
	discard_coords_append = discard_coords.append
	for line in cif:
		s = line.split()
		if isvert(s):
			discard_coords_append(s)

	start = np.asarray(map(float, discard_coords[0][2:5]))

	shortest_path_dict = nx.shortest_path(augTG)
	SN = sorted(augTG.nodes(), key = lambda x : int(re.sub('[A-Za-z]','',x)))
	sequential_paths = [(SN[i],SN[i+1]) for i in range(num_vertices) if i+1 < num_vertices]
	start = np.asarray(start)

	cnd = nx.get_node_attributes(augTG, 'cifname')
	coords = [[sequential_paths[0][0], cnd[sequential_paths[0][0]], start, [(e[2]['index'],e[2]['pd'],e[2]['cifname']) for e in augTG.edges(data=True) if sequential_paths[0][0] in e]]]
	coords_append = coords.append
	already_placed = [sequential_paths[0][0]]
	already_placed_append = already_placed.append

	for st in sequential_paths:
		path = shortest_path_dict[st[0]][st[1]]
		lp = len(path)
		traverse = [(path[i],path[i+1]) for i in range(lp) if i+1 < lp]
		for e0 in traverse:
			s,e = e0
			edict = augTG[s][e]
			key = [k for k in edict][0]
			ind = k[0]
			positive_direction = edict[key]['pd']

			if (s,e) == positive_direction:
				direction = 1
			elif (e,s) == positive_direction:
				direction = -1
			else:
				print 'Error in defining an edge traversal in omega_to_coords.py'
				print 'Exiting'
				sys.exit()
			start = start + direction * sc_omega_plus[ind - 1]

			if e0[1] not in already_placed:
				coords_append([e0[1], cnd[e0[1]], start, [(e[2]['index'], e[2]['pd'],e[2]['cifname']) for e in augTG.edges(data=True) if e0[1] in e]])
				already_placed_append(e0[1])
	
	norm_coords = []
	norm_coords_append = norm_coords.append
	for line in coords:
		vec = []
		vec_append = vec.append
		v = line[0]
		for dim in line[2]:
			if dim < 0.0:
				while dim < 0:
					dim = dim + 1.0
			elif dim >= 1.0:
				while dim >= 1.0:
					dim = dim - 1.0
			vec_append(dim)
		norm_coords_append([line[0],line[1],vec,line[3]])

	norm_coords = sorted(norm_coords, key = lambda x : int(re.sub('[A-Za-z]','',x[0])))

	path = os.path.join('templates', template)

	with open(path, 'r') as tcif:

		tcif = tcif.read()
		tcif = filter(None, tcif.split('\n'))

	if CHECK:

		q = num_vertices
		cpath = os.path.join('check_cifs', str(g) + '_check_scaled_' + template)

		with open(cpath, 'w') as check:
			for line in tcif:
				s = line.split()
				if not isvert(s):
					if '_cell_length_a' in line:
						check.write('_cell_length_a   ' + str(sc_a))
					elif '_cell_length_b' in line:
						check.write('_cell_length_b   ' + str(sc_b))
					elif '_cell_length_c' in line:
						check.write('_cell_length_c   ' + str(sc_c))
					elif '_cell_angle_alpha' in line:
						check.write('_cell_angle_alpha   ' + str(sc_alpha))
					elif '_cell_angle_beta' in line:
						check.write('_cell_angle_beta   ' + str(sc_beta))
					elif '_cell_angle_gamma' in line:
						check.write('_cell_angle_gamma   ' + str(sc_gamma))
					else:
						check.write(line)
					check.write('\n')
				else:
					for n in norm_coords:
						name = re.sub('[0-9]','',n[0])
						v = n[2]
						check.write('{:>5}{:>5}{:>20}{:>20}{:>20}{:>12}{:>8}{:>8}'.format(n[0],name,v[0],v[1],v[2],'0.00000','Uiso','1.00'))
						check.write('\n')
					break

	return norm_coords
