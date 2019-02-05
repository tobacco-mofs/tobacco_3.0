import numpy as np
import re
import networkx as nx
import sys

def remove_Fr(placed_all, bonds_all):
	
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

	return(new_placed_all, new_bonds_all)

	