import re
import networkx as nx

def remove_Fr(placed_all, bonds_all):

	G = nx.Graph()
	for n in placed_all:
		G.add_node(n[0])
	for l in bonds_all:
		G.add_edge(l[0],l[1])
	
	nconnections = []
	
	for l in placed_all:
	
		nbors = list(G.neighbors(l[0]))
		
		if 'Fr' in ''.join(nbors) and 'Fr' not in l[0]:
			count = len([nbor for nbor in nbors if 'Fr' in nbor])
			l[-3] = 'X' + re.sub('[^0-9]','',l[-3])
			nconnections.append([l[0],count])

	new_placed_all = []
	new_placed_all_append = new_placed_all.append

	new_bonds_all = []
	new_bonds_all_append = new_bonds_all.append

	for l in placed_all:
		if re.sub('[0-9]','',l[0]) != 'Fr':
			new_placed_all_append(l)

	for l in bonds_all:
		if re.sub('[0-9]','',l[0]) != 'Fr' and re.sub('[0-9]','',l[1]) != 'Fr':
			new_bonds_all_append(l)

	return(new_placed_all, new_bonds_all, nconnections)

	