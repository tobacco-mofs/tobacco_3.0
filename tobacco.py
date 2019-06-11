from reindex import apply_reindex
from ciftemplate2graph import ct2g
from vertex_edge_assign import vertex_assign, assign_node_vecs2edges
from cycle_cocyle import cycle_cocyle, Bstar_alpha
from bbcif_properties import cncalc
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
YOU_ARE_PATIENT = configuration.YOU_ARE_PATIENT

CHECK_NUMBER_OF_VERTICES_AND_EDGES = configuration.CHECK_NUMBER_OF_VERTICES_AND_EDGES
WRITE_CHECK_FILES = configuration.WRITE_CHECK_FILES
WRITE_CIF = configuration.WRITE_CIF

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
ORIENTATION_DEPENDENT_NODES = configuration.ORIENTATION_DEPENDENT_NODES # added as temporary fix to node orientation problem (effects only one node BB to date)

PLACE_EDGES_BETWEEN_CONNECTION_POINTS = configuration.PLACE_EDGES_BETWEEN_CONNECTION_POINTS
RECORD_CALLBACK = configuration.RECORD_CALLBACK
NET_2D = configuration.NET_2D
####### Global options #######

pi = np.pi

apply_reindex(CHARGES)

for d in ['templates', 'nodes', 'edges']:
	try:
		os.remove(os.path.join(d,'.DS_store'))
	except:
		pass

for template in os.listdir('templates'):
	print ''
	print '========================================================================================================='
	print '                                          template : ',template                                          
	print '========================================================================================================='
	print ''
	
	TG, unit_cell, TVT, TET, TNAME, a, b, c, ang_alpha, ang_beta, ang_gamma, max_le = ct2g(template)
	node_cns = [(cncalc(node, 'nodes', ONE_ATOM_NODE_CN), node) for node in os.listdir('nodes')]
	edge_type_key = dict((list(TET)[k],k) for k in xrange(len(TET)))

	print 'Number of vertices = ', len(TG.nodes())
	print 'Number of edges = ', len(TG.edges())
	print ''
	
	if PRINT:
		print '*****************************************************************'
		print '   There are', len(TG.nodes()), 'vertices in the voltage graph'
		print '*****************************************************************'
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
		print ''
		print '**************************************************************'
		print '   There are', len(TG.edges()), 'edges in the voltage graph'
		print '**************************************************************'
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

	vas = vertex_assign(TG, TVT, node_cns, unit_cell, ONE_ATOM_NODE_CN, USER_SPECIFIED_NODE_ASSIGNMENT, SYMMETRY_TOL)

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
		print ''
		print '************************************************************************'
		print '   Bstar (top) and alpha (bottom) for the barycentric embedding are : ' 
		print '************************************************************************'
		print ''
		for i in Bstar:
			print i
		print ''
		for i in alpha:
			print i
		print ''

	omega_plus = np.dot(np.linalg.inv(Bstar),alpha)
	num_vertices = len(TG.nodes())

	if CHECK_NUMBER_OF_VERTICES_AND_EDGES:
		defined = True
		try:
			rcsr_nv, rcsr_ne, rest = template.split('_')
		except ValueError:
			defined = False
			print 'Warning! edge and vertex numbers could not be checked.'
			print 'The template filename does not contain vertex and/or edge numbers.'
			print ''
		if defined:
			if int(rcsr_nv) != num_vertices:
				print 'Warning! the template does not have the correct number of vertices!'
			if int(rcsr_ne) != num_edges:
				print 'Warning! the template does not have the correct number of vertices!'

	if len([e for e in os.listdir('edges')]) == 0:
		print 'There are no edge cifs in the edge directory.'
		print 'Exiting'
		sys.exit()

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

		v_set = list(set([i[1] for i in va]))

		print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		print 'vertex assignment : ',v_set
		print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		print ''

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
	
			sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma,sc_covar,Bstar_inv,max_length,callbackresults,ncra,ncca = scale(all_SBU_coords,a,b,c,ang_alpha,ang_beta,ang_gamma,max_le,num_vertices,Bstar,alpha,num_edges,NET_2D,YOU_ARE_PATIENT,SCALING_ITERATIONS)
	
			print '*******************************************'
			print '   The scaled unit cell parameters are : ' 
			print '*******************************************'
			print ''
			print 'a : ', sc_a
			print 'b : ', sc_b
			print 'c : ', sc_c
			print 'alpha : ', sc_alpha
			print 'beta : ', sc_beta
			print 'gamma : ', sc_gamma
			print ''
	
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

				print callbackresults[-1][-1]

			if PLACE_EDGES_BETWEEN_CONNECTION_POINTS:
				placed_edges = adjust_edges(placed_edges, placed_nodes, sc_unit_cell)
	
			placed_all = placed_nodes + placed_edges
			bonds_all = node_bonds + edge_bonds
	
			if WRITE_CHECK_FILES:
				print 'writing check files...'
				print ''
				write_check_cif(template, placed_nodes, placed_edges, g, scaled_params, sc_unit_cell)
	
			if SINGLE_ATOM_NODE or NODE_TO_NODE:
				placed_all,bonds_all = remove_Fr(placed_all,bonds_all)
	
			print 'computing X-X bonds...'
			
			fixed_bonds, nbcount, bond_check = bond_connected_components(placed_all, bonds_all, sc_unit_cell, max_length, BOND_TOL, TRACE_BOND_MAKING, NODE_TO_NODE, EXPANSIVE_BOND_SEARCH, ONE_ATOM_NODE_CN)
	
			print ''
			print 'there were ' + str(nbcount) + ' X-X bonds formed'
			if bond_check:
				print 'bond check passed'
				bond_check_code = ''
			else:
				print 'bond check failed, attempting distance search bonding...'
				fixed_bonds, nbcount = distance_search_bond(placed_all, bonds_all, sc_unit_cell, 2.5, TRACE_BOND_MAKING)
				bond_check_code = '_BOND_CHECK'
				print 'there were', nbcount, 'X-X bonds formed'
			print ''
	
			if CHARGES:
				print 'fixing charges...'
				print ''
				fc_placed_all, netcharge, onetcharge, rcb = fix_charges(placed_all)
			else:
				fc_placed_all = placed_all

			fixed_bonds = fix_bond_sym(fixed_bonds, placed_all, sc_unit_cell)

			if CHARGES:
				print '*******************************************'
				print '            Charge information :           ' 
				print '*******************************************'
				print ''
				print 'old net charge                   : ', onetcharge
				print 'new net charge (after rescaling) : ', netcharge
				print 'rescaling magnitude              : ', rcb
				print ''
	
			vnames = '_'.join([v.split('.')[0] for v in v_set])
	
			if len(ea) <= 5:
				enames = '_'.join([e[0:-4] for e in ea])
			else:
				enames = str(len(ea)) + '_edges'
	
			cifname = template[0:-4] + '_' +  vnames + '_' + enames + bond_check_code + '.cif'
	
			if WRITE_CIF:
				print 'writing cif...'
				print ''
				write_cif(fc_placed_all, fixed_bonds, scaled_params, sc_unit_cell, cifname, CHARGES)

print 'Normal termination of Tobacco_3.0 after'
print '--- %s seconds ---' % (time.time() - start_time)
print ''
