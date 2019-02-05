import networkx as nx
import numpy as np
import re

def cycle_cocyle(TG):
	MST = nx.minimum_spanning_tree(TG)
	scaffold = nx.MultiGraph()

	used_keys = []
	for e in MST.edges(data=True):
		edict = e[2]
		lbl = edict['label']
		ind = edict['index']
		ke = (ind,lbl[0],lbl[1],lbl[2])
		scaffold.add_edge(e[0],e[1],key=ke)
		used_keys.append(ke)

	cycle_basis = []
	cycle_basis_append = cycle_basis.append
	nxfc = nx.find_cycle

	for e0 in TG.edges(data=True):
		edict = e0[2]
		lbl = edict['label']
		ind = edict['index']
		ke = (ind,lbl[0],lbl[1],lbl[2])
		scaffold.add_edge(e0[0],e0[1],key=ke)
		if ke not in used_keys:
			cycles = list(nxfc(scaffold))
			cy_list = [(i[0], i[1], i[2]) for i in cycles]
			if cy_list not in cycle_basis:
				cycle_basis_append(cy_list)
			scaffold.remove_edge(e0[0],e0[1],key=ke)

	node_out_edges = []
	node_out_edges_append = node_out_edges.append
	node_list = list(TG.nodes())
	
	for n in range(len(TG.nodes()) - 1):
		node = node_list[n]
		noe = [node]
		for e in TG.edges(data=True):
			if node == e[0] or node == e[1]:
				edict = e[2]
				lbl = edict['label']
				ke = (edict['index'],lbl[0],lbl[1],lbl[2])
				positive_direction = edict['pd']
				noe.append(positive_direction + ke)
		node_out_edges_append(noe)

	return cycle_basis, node_out_edges

def Bstar_alpha(CB, CO, augTG, num_edges):

	edge_keys = dict(((k[0],k[1],k[2]['index']),[]) for k in augTG.edges(data=True))
	for e in augTG.edges(keys=True, data=True):
		edge_keys[(e[0],e[1],e[3]['index'])] = e[2] 

	Bstar = []
	a = []
	Bstar_append = Bstar.append
	a_append = a.append
	q = 0
	for cycle in CB:
		q += 1
		cycle_vec = [0] * num_edges
		net_voltage = np.array([0,0,0])
		for edge in cycle:
			s,e,lv = edge
			ind = lv[0]
			voltage = np.asarray(lv[1:])
			try:
				key = edge_keys[(s,e,ind)]
			except:
				key = edge_keys[(e,s,ind)]
			positive_direction = augTG[s][e][key]['pd']
			if (s,e) == positive_direction:
				direction = 1
			elif (e,s) == positive_direction:
				direction = -1
			else:
				print 'Error in Bstar cycle vector construction, edge direction cannot be defined for:'
				print s,e

			cycle_vec[ind - 1] = direction
			net_voltage = net_voltage + (direction * voltage)

		Bstar_append(cycle_vec)
		a_append(net_voltage)

	for vertex in sorted(CO, key = lambda x: int(re.sub('[A-Za-z]','',x[0]))):
		cocycle_vec = [0] * num_edges
		v = vertex[0]
		ooa = [[i[2],(i[0],i[1]),np.array(i[3:])] for i in vertex[1:]]
		for out_edge in ooa:
			ke = (out_edge[0], out_edge[2][0], out_edge[2][1], out_edge[2][2])
			ind = out_edge[0]
			s,e = out_edge[1]
			positive_direction = augTG[s][e][ke]['pd'] 

			if '_a' not in s and '_a' not in e:
				if s == v:
					o = e
				else:
					o = s
				v_ind = int(re.sub('[A-Za-z]','',v))
				o_ind = int(re.sub('[A-Za-z]','',o))
				if v_ind < o_ind:
					cd = 1
				else:
					cd = -1
				if v == s:
					direction = 1
				else:
					direction = -1
				if direction != cd:
					print 'Warning! direction assignment for the co-cycle vector', s,e, 'may be incorrect.'
					print 'The direction assignment does not follow the low-index to high-index = positive convention'
			
			cocycle_vec[ind - 1] = direction

		Bstar_append(cocycle_vec)
		a_append(np.array([0,0,0]))

	if len(Bstar) != len(a):
		print 'Error in cycle_cocycle.py, the row ranks of Bstar and alpha do not match.'

	return np.asarray(Bstar), np.asarray(a)
