from reindex import apply_reindex
from ciftemplate2graph import ct2g
from vertex_edge_assign import vertex_assign, assign_node_vecs2edges
from cycle_cocyle import cycle_cocyle, Bstar_alpha
from bbcif_properties import cncalc, bbelems
from SBU_geometry import SBU_coords
from scale import scale
from scaled_embedding2coords import omega2coords
from place_bbs import scaled_node_and_edge_vectors, place_nodes, place_edges
from remove_net_charge import fix_charges
from remove_dummy_atoms import remove_Fr
from adjust_edges import adjust_edges
from write_cifs import write_check_cif, write_cif, bond_connected_components, distance_search_bond, fix_bond_sym
from scale_animation import scaling_callback_animation, write_scaling_callback_animation, animate_objective_minimization

import configuration
import glob
import os
import re
import networkx as nx
import numpy as np
import sys
import itertools
import time

start_time = time.time()

####### Global options #######
PRINT = configuration.PRINT
ONE_ATOM_NODE_CN = configuration.ONE_ATOM_NODE_CN
CONNECTION_SITE_BOND_LENGTH = configuration.CONNECTION_SITE_BOND_LENGTH
WRITE_CHECK_FILES = configuration.WRITE_CHECK_FILES
WRITE_CIF = configuration.WRITE_CIF
ALL_NODE_COMBINATIONS = configuration.ALL_NODE_COMBINATIONS
USER_SPECIFIED_NODE_ASSIGNMENT = configuration.USER_SPECIFIED_NODE_ASSIGNMENT
COMBINATORIAL_EDGE_ASSIGNMENT = configuration.COMBINATORIAL_EDGE_ASSIGNMENT
CHARGES = configuration.CHARGES
SCALING_ITERATIONS = configuration.SCALING_ITERATIONS
SYMMETRY_TOL = configuration.SYMMETRY_TOL
BOND_TOL = configuration.BOND_TOL
EXPANSIVE_BOND_SEARCH = configuration.EXPANSIVE_BOND_SEARCH
TRACE_BOND_MAKING = configuration.TRACE_BOND_MAKING
NODE_TO_NODE = configuration.NODE_TO_NODE
SINGLE_ATOM_NODE = configuration.SINGLE_ATOM_NODE
ORIENTATION_DEPENDENT_NODES = configuration.ORIENTATION_DEPENDENT_NODES
PLACE_EDGES_BETWEEN_CONNECTION_POINTS = configuration.PLACE_EDGES_BETWEEN_CONNECTION_POINTS
RECORD_CALLBACK = configuration.RECORD_CALLBACK
OUTPUT_SCALING_DATA = configuration.OUTPUT_SCALING_DATA
FIX_UC = configuration.FIX_UC
PRE_SCALE = configuration.PRE_SCALE
SCALING_CONVERGENCE_TOLERANCE = configuration.SCALING_CONVERGENCE_TOLERANCE
SCALING_STEP_SIZE = configuration.SCALING_STEP_SIZE
SINGLE_METAL_MOFS_ONLY = configuration.SINGLE_METAL_MOFS_ONLY
####### Global options #######

pi = np.pi

vname_dict = {'V':1  ,'Er':2 ,'Ti':3 ,'Ce':4 ,'S':5  ,
			  'H':6  ,'He':7 ,'Li':8 ,'Be':9 ,'B':10 ,
			  'C':11 ,'N':12 ,'O':13 ,'F':14 ,'Ne':15,
			  'Na':16,'Mg':17,'Al':18,'Si':19,'P':20 ,
			  'Cl':21,'Ar':22,'K':23 ,'Ca':24,'Sc':24,
			  'Cr':26,'Mn':27,'Fe':28,'Co':29,'Ni':30 }

metal_elements = ['Cu', 'Cr', 'Zn', 'Zr', 'Fe', 'Al']

#apply_reindex(CHARGES)

for d in ['templates', 'nodes', 'edges']:
	try:
		os.remove(os.path.join(d,'.DS_store'))
	except:
		pass

if OUTPUT_SCALING_DATA:
	scd = open('scaling_data.txt', 'w')
	scd.write('ab_ratio_i ab_ratio_f ac_ratio_i ac_ratio_f bc_ratio_i bc_ratio_f alpha_i alpha_f beta_i beta_f gamma_i gamma_f average_covec final_obj\n')

for template in sorted(os.listdir('templates')):
	print ''
	print '========================================================================================================='
	print 'template :',template                                          
	print '========================================================================================================='
	print ''
	
	TG, unit_cell, TVT, TET, TNAME, a, b, c, ang_alpha, ang_beta, ang_gamma, max_le = ct2g(template)
	node_cns = [(cncalc(node, 'nodes', ONE_ATOM_NODE_CN), node) for node in os.listdir('nodes')]
	edge_type_key = dict((list(TET)[k],k) for k in xrange(len(TET)))

	print 'Number of vertices = ', len(TG.nodes())
	print 'Number of edges = ', len(TG.edges())
	print ''
	
	if PRINT:
		print 'There are', len(TG.nodes()), 'vertices in the voltage graph:'
		print ''
		q = 0
		for node in TG.nodes():
			q += 1
			print q,':',node
			node_dict = TG.node[node]
			print 'type : ', node_dict['type']
			print 'cartesian coords : ', node_dict['ccoords']
			print 'fractional coords : ', node_dict['fcoords']
			print 'degree : ', node_dict['cn'][0]
			print ''
		print 'There are', len(TG.edges()), 'edges in the voltage graph:'
		print ''
		for edge in TG.edges(data=True,keys=True):
			edge_dict = edge[3]
			ind = edge[2]
			print ind,':',edge[0],edge[1]
			print 'length : ',edge_dict['length']
			print 'type : ',edge_dict['type']
			print 'label : ',edge_dict['label']
			print 'positive direction :',edge_dict['pd']
			print 'cartesian coords : ',edge_dict['ccoords']
			print 'fractional coords : ',edge_dict['fcoords']
			print ''

	vas = vertex_assign(TG, TVT, node_cns, unit_cell, ONE_ATOM_NODE_CN, USER_SPECIFIED_NODE_ASSIGNMENT, SYMMETRY_TOL, ALL_NODE_COMBINATIONS)

	num_edges = len(TG.edges())
	CB,CO = cycle_cocyle(TG)
	
	for va in vas:
		if len(va) == 0:
			print 'At least one vertex does not have a building block with the correct number of connection sites.'
			print 'Moving to the next template.'
			print ''
			continue

	if len(CB) != (len(TG.edges()) - len(TG.nodes()) + 1):
		print 'The cycle basis is incorrect.'
		print 'The number of cycles in the cycle basis does not equal the rank of the cycle space.'
		print 'Moving to the next tempate.'
		continue

	Bstar, alpha  = Bstar_alpha(CB,CO,TG,num_edges)

	if PRINT:
		print 'star (top) and alpha (bottom) for the barycentric embedding are: ' 
		print ''
		for i in Bstar:
			print i
		print ''
		for i in alpha:
			print i
		print ''

	omega_plus = np.dot(np.linalg.inv(Bstar),alpha)
	num_vertices = len(TG.nodes())

	if COMBINATORIAL_EDGE_ASSIGNMENT:
		eas = list(itertools.product([e for e in os.listdir('edges')], repeat = len(TET)))
	else:
		edge_files = sorted([e for e in os.listdir('edges')])
		eas = []
		i = 0
		while len(eas) < len(TET):
			eas.append(edge_files[i])
			i += 1
			if i == len(edge_files):
				i = 0
		eas = [eas]

	g = 0
	for va in vas:

		node_elems = [bbelems(i[1], 'nodes') for i in va]
		metals = [[i for i in j if i in metal_elements] for j in node_elems]
		metals = set([i for j in metals for i in j])

		v_set = [('v' + str(vname_dict[re.sub('[0-9]','',i[0])]), i[1]) for i in va]
		v_set = sorted(list(set(v_set)), key=lambda x: x[0])
		v_set = [v[0] + '-' + v[1] for v in v_set]

		print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		print 'vertex assignment : ',v_set
		print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		print ''

		if len(metals) != 1 and SINGLE_METAL_MOFS_ONLY:
			print v_set, 'contains no metals or multiple metal elements, no cif will be written', metals
			print ''
			continue

		for v in va:
			for n in TG.nodes(data=True):
				if v[0] == n[0]:
					n[1]['cifname'] = v[1]

		for ea in eas:

			g += 1

			print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
			print 'edge assignment : ',ea
			print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
			print ''

			type_assign = dict((k,[]) for k in TET)
			for k,m in zip(TET,ea):
				type_assign[k] = m
	
			for e in TG.edges(data=True):
				ty = e[2]['type']
				for k in type_assign:
					if ty == k or (ty[1],ty[0]) == k:
						e[2]['cifname'] = type_assign[k]
	
			ea_dict = assign_node_vecs2edges(TG, unit_cell, SYMMETRY_TOL)
			all_SBU_coords = SBU_coords(TG, ea_dict, CONNECTION_SITE_BOND_LENGTH)
			sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma, sc_covar, Bstar_inv, max_length, callbackresults, ncra, ncca, scaling_data = scale(all_SBU_coords,a,b,c,ang_alpha,ang_beta,ang_gamma,max_le,num_vertices,Bstar,alpha,num_edges,FIX_UC,SCALING_ITERATIONS,PRE_SCALE,SCALING_CONVERGENCE_TOLERANCE,SCALING_STEP_SIZE)
	
			print '*******************************************'
			print 'The scaled unit cell parameters are : ' 
			print '*******************************************'
			print 'a    :', np.round(sc_a, 5)
			print 'b    :', np.round(sc_b, 5)
			print 'c    :', np.round(sc_c, 5)
			print 'alpha:', np.round(sc_alpha, 5)
			print 'beta :', np.round(sc_beta, 5)
			print 'gamma:', np.round(sc_gamma, 5)
			print ''

			for sc, name in zip((sc_a, sc_b, sc_c), ('a', 'b', 'c')):
				cflag = False
				if sc < 1.0:
					print 'unit cell parameter', name, 'has collapsed during scaling!'
					print 'try re-running with', name, 'fixed, with a larger value for PRE_SCALE, or with a higher SCALING_CONVERGENCE_TOLERANCE'
					print 'no cif will be written'
					cflag = True

			if cflag:
				continue

			scaled_params = [sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma]
		
			sc_Alpha = np.r_[alpha[0:num_edges-num_vertices+1,:], sc_covar]
			sc_omega_plus = np.dot(Bstar_inv, sc_Alpha)
		
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
			
			scaled_coords = omega2coords(TG, sc_omega_plus, (sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma), num_vertices, template, g, WRITE_CHECK_FILES)
			nvecs,evecs = scaled_node_and_edge_vectors(scaled_coords, sc_omega_plus, sc_unit_cell, ea_dict)
			placed_nodes, node_bonds = place_nodes(nvecs, CHARGES, ORIENTATION_DEPENDENT_NODES)
			placed_edges, edge_bonds = place_edges(evecs, CHARGES, len(placed_nodes))

			if RECORD_CALLBACK:

				vnames = '_'.join([v.split('.')[0] for v in v_set])

				if len(ea) <= 5:
					enames = '_'.join([e[0:-4] for e in ea])
				else:
					enames = str(len(ea)) + '_edges'

				prefix = template[0:-4] + '_' +  vnames + '_' + enames

				frames = scaling_callback_animation(callbackresults, alpha, Bstar_inv, ncra, ncca, num_vertices, num_edges, TG, template, g, False)
				write_scaling_callback_animation(frames, prefix)
				animate_objective_minimization(callbackresults, prefix)

			if PLACE_EDGES_BETWEEN_CONNECTION_POINTS:
				placed_edges = adjust_edges(placed_edges, placed_nodes, sc_unit_cell)
	
			placed_all = placed_nodes + placed_edges
			bonds_all = node_bonds + edge_bonds
	
			if WRITE_CHECK_FILES:
				write_check_cif(template, placed_nodes, placed_edges, g, scaled_params, sc_unit_cell)
	
			if SINGLE_ATOM_NODE or NODE_TO_NODE:
				placed_all,bonds_all = remove_Fr(placed_all,bonds_all)
			
			print 'computing X-X bonds...'
			print ''
			print '*******************************************'
			print 'Bond formation : ' 
			print '*******************************************'
			
			fixed_bonds, nbcount, bond_check = bond_connected_components(placed_all, bonds_all, sc_unit_cell, max_length, BOND_TOL, TRACE_BOND_MAKING, NODE_TO_NODE, EXPANSIVE_BOND_SEARCH, ONE_ATOM_NODE_CN)

			print 'there were ' + str(nbcount) + ' X-X bonds formed'
			if bond_check:
				print 'bond check passed'
				bond_check_code = ''
			else:
				continue
				print 'bond check failed, attempting distance search bonding...'
				fixed_bonds, nbcount = distance_search_bond(placed_all, bonds_all, sc_unit_cell, 2.5, TRACE_BOND_MAKING)
				bond_check_code = '_BOND_CHECK'
				print 'there were', nbcount, 'X-X bonds formed'
			print ''
	
			if CHARGES:
				fc_placed_all, netcharge, onetcharge, rcb = fix_charges(placed_all)
			else:
				fc_placed_all = placed_all

			fixed_bonds = fix_bond_sym(fixed_bonds, placed_all, sc_unit_cell)

			if CHARGES:
				print '*******************************************'
				print 'Charge information :                       ' 
				print '*******************************************'
				print 'old net charge                  :', np.round(onetcharge, 5)
				print 'new net charge (after rescaling):', np.round(netcharge, 5)
				print 'rescaling magnitude             :', np.round(rcb, 5)
				print ''
	
			vnames = '_'.join([v.split('.')[0] for v in v_set])
	
			if len(ea) <= 5:
				enames = []
				for e in [e[0:-4] for e in ea]:
					if e not in enames:
						enames.append(e)
				enames = '_'.join(enames)

			else:
				enames = str(len(ea)) + '_edges'
	
			cifname = template[0:-4] + '_' +  vnames + '_' + enames + bond_check_code + '.cif'

			if OUTPUT_SCALING_DATA:
				line = [cifname.split('.')[0]] + [i for j in scaling_data for i in j]
				format_string = ' '.join(['{}' for i in line])
				scd.write(format_string.format(*line))
				scd.write('\n')
	
			if WRITE_CIF:
				print 'writing cif...'
				print ''
				write_cif(fc_placed_all, fixed_bonds, scaled_params, sc_unit_cell, cifname, CHARGES)

if OUTPUT_SCALING_DATA:
	scd.close()
	
print 'Normal termination of Tobacco_3.0 after'
print '--- %s seconds ---' % (time.time() - start_time)
print ''
