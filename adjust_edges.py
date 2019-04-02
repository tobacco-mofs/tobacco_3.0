import numpy as np
import re 

from place_bbs import superimpose

def PBC3DF_sym(vec1, vec2):

	dX,dY,dZ = vec1 - vec2
			
	if dX > 0.5:
		s1 = 1.0
		ndX = dX - 1.0
	elif dX < -0.5:
		s1 = -1.0
		ndX = dX + 1.0
	else:
		s1 = 0.0
		ndX = dX
				
	if dY > 0.5:
		s2 = 1.0
		ndY = dY - 1.0
	elif dY < -0.5:
		s2 = -1.0
		ndY = dY + 1.0
	else:
		s2 = 0.0
		ndY = dY
	
	if dZ > 0.5:
		s3 = 1.0
		ndZ = dZ - 1.0
	elif dZ < -0.5:
		s3 = -1.0
		ndZ = dZ + 1.0
	else:
		s3 = 0.0
		ndZ = dZ

	sym = np.array([s1,s2,s3])

	return np.array([ndX,ndY,ndZ]), sym

def adjust_edges(placed_edges, placed_nodes, sc_unit_cell):

	adjusted_placed_edges = []
	adjusted_placed_edges_extend = adjusted_placed_edges.extend

	placed_edges = np.asarray(placed_edges)
	edge_labels = set(map(int, placed_edges[:,-1]))

	edge_dict = dict((k,[]) for k in edge_labels)

	node_connection_points = [map(float,i[1:4]) for i in placed_nodes if re.sub('[0-9]','',i[5]) == 'X']

	for edge in placed_edges:
		ty = int(edge[-1])
		edge_dict[ty].append(edge)

	for k in edge_dict:

		edge = np.asarray(edge_dict[k])
		elems = edge[:,0]
		evecs = [map(float,i) for i in edge[:,1:4]]
		charges = edge[:,4]
		cp = edge[:,5]
		ty = edge[:,6]

		xvecs = [map(float,i) for (i,j) in zip(evecs,cp) if re.sub('[0-9]','',j) == 'X']
		relavent_node_xvecs = []
		relavent_node_xvecs_append = relavent_node_xvecs.append

		for ex in xvecs:

			min_dist = (1e6, [], 0)

			f_ex = np.dot(np.linalg.inv(sc_unit_cell), ex)
			for i in range(len(node_connection_points)):
				nx = node_connection_points[i]
				f_nx = np.dot(np.linalg.inv(sc_unit_cell), nx)

				fdist_vec,sym = PBC3DF_sym(f_ex, f_nx) 
				cdist = np.linalg.norm(np.dot(sc_unit_cell, fdist_vec))

				if cdist < min_dist[0]:
					min_dist = (cdist, np.dot(sc_unit_cell,f_nx + sym), i)

			node_connection_points.pop(min_dist[2])
			relavent_node_xvecs_append(min_dist[1])

		ecom = np.average(xvecs, axis=0)
		rnxcom = np.average(relavent_node_xvecs, axis=0)

		evecs = np.asarray(evecs - ecom)
		xvecs = np.asarray(xvecs - ecom)
		relavent_node_xvecs = np.asarray(relavent_node_xvecs)

		trans = rnxcom
		min_dist,rot,tran = superimpose(xvecs,relavent_node_xvecs)
		adjusted_evecs = np.dot(evecs,rot) + trans
		adjusted_edge = np.column_stack((elems,adjusted_evecs,charges,cp,ty))
		adjusted_placed_edges_extend(adjusted_edge)

	return(adjusted_placed_edges)		
