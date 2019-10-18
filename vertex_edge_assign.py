import os
import itertools
import numpy as np
from bbcif_properties import X_vecs
from place_bbs import superimpose
from ciftemplate2graph import node_vecs

def vertex_assign(TG, TVT, node_cns, unit_cell, cn1, USNA, SYM_TOL, ALL_NODE_COMBINATIONS):

	node_dict = dict((k,[]) for k in TVT)

	for node in node_cns:
		for k in TVT:
			if node[0] == k[0]:
				node_dict[k].append(node[1])

	if USNA:

		va = []
		va_append = va.append

		choice_dict = dict((k,'') for k in TVT)
		if not os.path.isfile('vertex_assignment.txt'):
			for k in node_dict:
				cn,name = k
				print ''
				print '???????????????????????????????????'
				print 'select building block for:', name , '(CN=' + str(cn) + ')'
				for c in range(len(node_dict[k])):
					print c, node_dict[k][c]
				cif_index = int(raw_input('enter the index of the desired cif: \n'))
				choice_dict[k] = node_dict[k][cif_index]
				print '???????????????????????????????????'
				print ''
		else:
			with open('vertex_assignment.txt','r') as va_key:
				va_key = va_key.read()
				va_key = va_key.split('\n')
				choices = [(l.split()[0],l.split()[1]) for l in va_key if len(l.split())==2]
			for k in node_dict:
				for c in choices:
					if c[0] == k[1] and c[1] in node_dict[k]:
						choice_dict[k] = c[1]
						break
					else:
						continue

		for k in choice_dict:

			if len(choice_dict[k]) == 0:
				print 'Error in vertex_edge_assign.py:'
				print 'Node type', k[0], 'has not assigned cif.'
				print 'Exiting'
				sys.exit()

			for n in TG.nodes(data=True):
				name,ndict = n
				if ndict['type'] == k[1]:
					va_append((name, choice_dict[k]))

		va = [va]

	else:

		print '*****************************************************************'
		print 'RMSD of the compatible node BBs with assigned vertices:          '
		print '*****************************************************************'
		print ''
		
		RMSDs = []
		RMSDs_append = RMSDs.append
		sym_assign = []
		sym_assign_append = sym_assign.append

		for k in node_dict:

			print 'vertex', k[1], '('+str(k[0]) + ' connected)'

			matched = 0
			unmatched = 0
			
			if len(node_dict[k]) == 0:
				continue
			coord_num = k[0]
	
			for n in TG.nodes(data=True):
				name,ndict = n
				distances = []
				distances_append = distances.append

				if ndict['type'] == k[1]:
					for cif in node_dict[k]:

						nvec = np.array([v/np.linalg.norm(v) for v in node_vecs(name, TG, unit_cell, False)])
						bbxvec = np.array([v/np.linalg.norm(v) for v in X_vecs(cif, 'nodes', False)])
						rmsd,rot,tran = superimpose(bbxvec,nvec)
						aff_b = np.dot(bbxvec,rot) + tran
						distances_append((rmsd,cif))

					for d in distances:
						disp,cif = d
						if d[0] < SYM_TOL[coord_num]:
							matched += 1
							matches = '(within tolerance)'
						else:
							unmatched += 1
							matches = '(outside tolerance)'
						print '    ', cif, 'deviation =', np.round(disp,5), matches

					for d in distances:
						if d[0] < SYM_TOL[coord_num]:
							sym_assign_append((k[1],d[1]))
					break
			print '*', matched, 'compatible building blocks out of', len(node_dict[k]), 'available for node', k[1], '*'
		print ''
		
		rearrange = dict((k[1],[]) for k in TVT)
		for a in sym_assign:
			rearrange[a[0]].append((a[0],a[1]))

		va_uncomb = [rearrange[a] for a in rearrange]

		va = []
		va_append = va.append
		used = []
		used_append = used.append
		for l in itertools.product(*va_uncomb):

			cifs = sorted(tuple([c[1] for c in l]))
			if cifs in used and not ALL_NODE_COMBINATIONS:
				continue

			choice_dict = dict((i[0],i[1]) for i in l)
			va_temp = []
			va_temp_append = va_temp.append
			for n in TG.nodes(data=True):
				name,ndict = n
				va_temp_append((name, choice_dict[ndict['type']]))
			va_append(va_temp)
			used_append(cifs)
					
	return va

def assign_node_vecs2edges(TG, unit_cell, SYM_TOL):
	edge_assign_dict = dict((k,{}) for k in TG.nodes())

	for n in TG.nodes(data=True):
		name,ndict = n
		cif = ndict['cifname']
		bbx = X_vecs(cif, 'nodes', True)
		nod = node_vecs(n[0], TG, unit_cell, True)

		bbxlabels = np.array([l[0] for l in bbx])
		nodlabels = np.array([l[0] for l in nod])

		bbxvec = np.array([l[1] for l in bbx])
		
		ll = 0
		for v in bbxvec:
			mag = np.linalg.norm(v - np.average(bbxvec, axis = 0))
			if mag > ll:
				ll = mag

		nodvec = np.array([mag * (l[1]/np.linalg.norm(l[1])) for l in nod])

		rmsd,rot,tran = superimpose(bbxvec, nodvec)

		aff_b = np.dot(bbxvec,rot) + tran

		laff_b = np.c_[bbxlabels,aff_b]
		lnodvec = np.c_[nodlabels,nodvec]
		
		asd = []
		asd_append = asd.append
		for v1 in laff_b:
			smallest_dist = (1.0E6, 'foo', 'bar')
			v1vec = map(float, v1[1:])
			mag = np.linalg.norm(v1vec)
			for v2 in lnodvec:
				dist = np.linalg.norm(v1vec - v2[1:])
				if dist < smallest_dist[0]:
					smallest_dist = (dist, v1[0], int(v2[0]), mag, v1vec)

			asd_append(smallest_dist)
			
		elad = dict((k[2],(k[1], k[3], k[4])) for k in asd)
		edge_assign_dict[name] = elad
					
	return edge_assign_dict
	
