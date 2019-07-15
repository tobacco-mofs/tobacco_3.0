import re
import os
import sys
import numpy as np
import networkx as nx

vertices = ('V','Er','Ti','Ce','S',
			'H','He','Li','Be','B',
			'C','N','O','F','Ne',
			'Na','Mg','Al','Si','P',
			'Cl','Ar','K','Ca','Sc',
			'Cr','Mn','Fe','Co','Ni')
pi = np.pi

def isfloat(value):
	
	try:
		float(value)
		return True
	except ValueError:
		return False

def isvert(line):

	if len(line) >= 5:
		if re.sub('[^a-zA-Z]','', line[0]) in vertices and line[1] in vertices:
			if (isfloat(line[2]) and isfloat(line[3]) and isfloat(line[4])):
				return True
			else:
				return False
		elif re.sub('[^a-zA-Z]','', line[0]) not in vertices or line[1] not in vertices:
			if (isfloat(line[2]) and isfloat(line[3]) and isfloat(line[4])):
				print 'Error in ciftemplate2graph, invalid vertex name'
				print line
				print 'Exiting'
				sys.exit()
	else:
		return False
	
def isedge(line):

	if len(line) >= 5:
		if (re.sub('[^a-zA-Z]','', line[0]) in vertices and
		re.sub('[^a-zA-Z]','', line[1]) in vertices and 
		isfloat(line[2]) and not isfloat(line[3]) and
		not isfloat(line[4])):
			return True
		else:
			return False

def PBC3DF(c1, c2):

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

def ct2g(cifname):

	path = os.path.join('templates', cifname)

	with open(path, 'r') as template:
		template = template.read()
		template = filter(None, template.split('\n'))

	DG = nx.MultiDiGraph()
	G = nx.MultiGraph()

	edge_exist = False
	for line in template:
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

	nc = 0
	ne = 0

	types = []
	e_types = []
	G_forward_edges = []
	aae = []

	G_forward_edges_append = G_forward_edges.append
	types_append = types.append
	e_types_append = e_types.append
	aae_append = aae.append

	max_le = 1.0e6

	for line in template:
		s = line.split()
		ty = re.sub('[^a-zA-Z]','',s[0])
		types_append(ty)
		if isvert(s):
			nc += 1
			f_nvec = np.asarray(map(float, s[2:5]))
			c_nvec = np.dot(unit_cell, f_nvec)
			G.add_node(s[0], type=ty, index=nc, ccoords=c_nvec, fcoords=f_nvec, cn=[], cifname=[])
		if isedge(s):
			edge_exist = True

			l = sorted([re.sub('[^a-zA-Z]','',s[0]),re.sub('[^a-zA-Z]','',s[1])])
			e_types_append((l[0],l[1]))

			if '_' in s[3]:
				lbl = np.asarray(map(int, s[3].split('_')[1])) - 5
			elif s[3] == '.':
				lbl = np.array([0,0,0])
			else:
				print 'Error in ciftemplate2graph, there are unrecognized bond translational symmetries in', cifname
				print 'Exiting'
				sys.exit()
			nlbl = -1*lbl

			if (
			   (s[0],s[1],lbl[0],lbl[1],lbl[2]) not in aae and
			   (s[1],s[0],lbl[0],lbl[1],lbl[2]) not in aae and
			   (s[0],s[1],nlbl[0],nlbl[1],nlbl[2]) not in aae and
			   (s[1],s[0],nlbl[0],nlbl[1],nlbl[2]) not in aae
			   ):

				ne += 1
				aae_append((s[0],s[1],lbl[0],lbl[1],lbl[2]))

				v1 = G.node[s[0]]['fcoords']
				v2 = G.node[s[1]]['fcoords'] + lbl

				ef_coords = np.average(np.array([v1, v2]), axis=0)
				ec_coords = np.dot(unit_cell, ef_coords)

				cdist = np.linalg.norm(np.dot(unit_cell, v1 - v2))

				le = float(s[2])
				if cdist < max_le:
					max_le = cdist

				ke = (ne,lbl[0],lbl[1],lbl[2])
				
				G.add_edge(s[0],s[1], key=(ne,lbl[0],lbl[1],lbl[2]), label=lbl , length=le, fcoords=ef_coords, ccoords=ec_coords, index=ne, pd=(s[0],s[1]))

	e_types = set(e_types)
	for edge in G.edges(data=True):
		ni = (re.sub('[^a-zA-Z]','',edge[0]),re.sub('[^a-zA-Z]','',edge[1]))
		edict = edge[2]
		for et in e_types:
			if ni == et or ni[::-1] == et:
				edict['type'] = et

	cns = []
	cns_append = cns.append
	for node in G.nodes():
		ndict = G.node[node]
		cn = G.degree(node)
		ndict['cn'].append(cn)
		cns_append((cn, ndict['type']))
	cns = set(cns)

	if not edge_exist:
		print 'Error in ciftemplate2graph, no edges are given in the template:', cifname
		print 'Exiting'
		sys.exit()
		
	return (G, unit_cell, cns, e_types, cifname, a, b, c, alpha, beta, gamma, max_le)

def node_vecs(node, G, unit_cell, label):

	edge_coords = []
	edge_coords_append = edge_coords.append
	for e in G.edges(data=True):

		edict = e[2]
		positive_direction = edict['pd']
		lbl = edict['label']
		ind = edict['index']

		if node in e[0:2]:

			if node == positive_direction[0]:
				vec = edict['ccoords']
			else:
				vec = np.dot(unit_cell, -1 * lbl + edict['fcoords'])
			if label:
				edge_coords_append([ind, vec])
			else:
				edge_coords_append(vec)

	if label:
		ec_com = np.average(np.asarray([v[1] for v in edge_coords]), axis=0)
	else:
		ec_com = np.average(edge_coords, axis = 0)
	if label:
		shifted_edge_coords = [[v[0], v[1] - ec_com] for v in edge_coords]
	else:
		shifted_edge_coords = [vec - ec_com for vec in edge_coords]

	return shifted_edge_coords

def edge_vecs(edge, G, unit_cell):

	for e in G.edges(data=True):
		edict = e[2]
		if edict['index'] == edge:
			s,e = e[0:2]
			ccoords = edict['ccoords']
			v1 = G.node[s]['ccoords']
			v2 = G.node[e]['ccoords']

			return [v1 - ccoords, v2 - ccoords]