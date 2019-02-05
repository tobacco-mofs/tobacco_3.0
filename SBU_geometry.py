from bbcif_properties import X_vecs, calc_edge_len
import numpy as np

def SBU_coords(augTG, ea_dict, csbl):

	SBU_coords = []
	SBU_coords_append = SBU_coords.append
	for node in augTG.nodes(data=True):
		vertex = node[0]
		ndict = node[1]
		ncif = ndict['cifname']
		
		xvecs = []
		xvecs_append = xvecs.append
		for e in augTG.edges(data=True):
			edict = e[2]
			if vertex in e:

				ecif = edict['cifname']
				positive_direction = edict['pd']
				ind = edict['index']
				length = calc_edge_len(ecif,'edges')

				if vertex == positive_direction[0]:
					direction = 1
					ov = positive_direction[1]
				else:
					direction = -1
					ov = positive_direction[0]

				xvecname,dx_v,xvec = ea_dict[vertex][ind]
				dx_ov = ea_dict[ov][ind][1]

				if length < 0.1:
					total_length = dx_v + dx_ov + csbl
				else:
					total_length = dx_v + dx_ov + length + 2*csbl
				
				svec = (xvec/np.linalg.norm(xvec)) * total_length * direction
				xvecs_append([ind, svec])

		SBU_coords_append((vertex, xvecs))

	return SBU_coords
