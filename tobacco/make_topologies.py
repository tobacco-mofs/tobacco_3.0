import pymatgen as pm
from ase.geometry.cell import cellpar_to_cell
import os
from datetime import datetime
import numpy as np
import warnings

#SETTINGS BY USER
tol = 1E-2 #tolerance for distances
scale = 10 #scale lattice constants by this factor
cgd_filename = 'RCSRnets-2019-06-01.cgd' #http://rcsr.anu.edu.au/systre
consider_2D = False #consider 2D topologies (disabled by default)

#Internal settings
vnames = [
    'V','Er','Ti','Ce','S',
    'H','He','Li','Be','B',
    'C','N','O','F','Ne',
    'Na','Mg','Al','Si','P',
    'Cl','Ar','K','Ca','Sc',
    'Cr','Mn','Fe','Co','Ni'] #names of vertices
edge_center_name = 'Lr' #placeholder edge name

if edge_center_name in vnames:
    raise ValueError('Edge center name must not be in vnames',edge_center_name)

#List of 2D topologies where a is the dummy axis
dummya_list = ['cpr','cqx','sdd','sdf','sdh','sdi','sdo','sdv','sdz','tdv','tdz']

#List of 2D topologies where b is the dummy axis
dummyb_list = ['cqe','cqv','dhb','krv','krvd','krw','krwd','sdc','sdm','sdp','sdq','sdw','sdy','tdr','tdw','tdx','tdy']

#initialize lists
topologies_all = [] #all topologies
groups_all = [] #all spacegroups
cellpars_all = [] #all [a,b,c,alpha,beta,gamma]
vertices_all = [] #all [x,y,z] fractional positions of vertices
edges_center_all = [] #all [x,y,z] fractional positions of edge centers
edges_head_all = [] #all [x,y,z] fractional positions of edge heads
edges_tail_all = [] #all [x,y,z] fractional positions of edge tails
cn_all = [] #all vertex coordination numbers, coded as dictionaries
is_threedim_all = [] #all booleans for if 3D

#Make sure .cgd file is present
if not os.path.exists(cgd_filename):
    raise ValueError('Missing RCSR .cgd data file', cgd_filename)

#Forbidden names on Windows
forbidden_names = ['con','prn','aux','nul']

#Read info from .cgd file
with open(cgd_filename,'r') as r:
    for line in r:
        line = line.strip()

        #Initialize values for new topology
        if 'crystal' in line.lower():
            three_dim = True
            bad = False
            vertices = []
            edges_center = []
            edges_head = []
            edges_tail = []
            cn = {}
            vertices_count = 0

        #Get the topology name
        elif 'name' in line.lower():
            topology_val = line.lower().split('name')[-1].replace('*','_star').replace('-','').strip()
            if topology_val in forbidden_names:
                topology_val += '0'

        #Get the spacegroup
        elif 'group' in line.lower():

            #Do not alter capitalization of spacegroups
            group_val = line.split('GROUP')[-1].split('group')[-1].strip()

            #Use updated group names
            if group_val == 'Cmca':
                group_val = 'Cmce'

            #2D-->3D group names
            elif group_val == 'p4gm':
                group_val = 'P4/mbm'
            elif group_val == 'p4mm':
                group_val = 'P4/mmm'
            elif group_val == 'p6mm':
                group_val = 'P6/mmm'
            elif group_val == 'p4mm':
                group_val = 'P4/mmm'
            elif group_val == 'p6':
                group_val = 'P6/m'
            elif group_val == 'c2mm':
                group_val = 'Cmmm'
            elif group_val == 'p31m':
                group_val = 'P-62m'
            elif group_val == 'p2mg':
                group_val = 'Pmma'
            elif group_val == 'p3m1':
                group_val = 'P-6m2'
            elif group_val == 'p2gg':
                group_val = 'Pbam'
            elif group_val == 'p2mm':
                group_val = 'Pmmm'
            elif group_val == 'cm':
                group_val = 'Amm2'
            elif group_val == 'pg':
                group_val = 'Pmc21'
            elif group_val == 'p1':
                group_val = 'Pm'
            elif group_val == 'p2':
                group_val = 'P2/m'
            elif group_val == 'p2gg':
                group_val = 'Pbam'

        #Get the lattice constants (with scale*(a,b,c))
        elif 'cell' in line.lower():
            cell_val = line.lower().split('cell')[-1]
            cell_val = [float(i) for i in cell_val.split()]

            #Add dummy dimensions for 2D topologies
            if len(cell_val) == 3:
                three_dim = False
                if topology_val in dummya_list:
                    cell_alpha_temp = cell_val[2]
                    del cell_val[2]
                    cell_val.insert(0,2.0)
                    cell_val.extend([cell_alpha_temp,90.0,90.0])                    
                elif topology_val in dummyb_list:
                    cell_beta_temp = cell_val[2]
                    del cell_val[2]
                    cell_val.insert(1,2.0)
                    cell_val.extend([90.0,cell_beta_temp,90.0])    
                else:
                    cell_gamma_temp = cell_val[2]
                    del cell_val[2]
                    cell_val.append(2.0)
                    cell_val.extend([90.0,90.0,cell_gamma_temp])

            if len(cell_val) != 6:
                bad = True
                continue

            cell_val[0] = cell_val[0]*scale
            cell_val[1] = cell_val[1]*scale
            cell_val[2] = cell_val[2]*scale

        #Get the vertices (make sure it's 3D) and get CNs
        elif 'node' in line.lower() or 'atom' in line.lower():
            vert_val = line.lower().split('node')[-1].split('atom')[-1].strip()
            vert_val = [i for i in vert_val.split()]

            #Add 0.0 dummy coordinate for 2D topologies
            if len(vert_val) == 4:
                if topology_val in dummya_list:
                    vert_val.insert(2,'0.0')
                elif topology_val in dummyb_list:
                    vert_val.insert(3,'0.0')
                else:
                    vert_val.append('0.0')

            #Make sure there is a coordination number and [x,y,z]
            if len(vert_val) != 5:
                bad = True
                continue

            vertices.append([float(vert_val[2]),float(vert_val[3]),float(vert_val[4])])
            vertices_count += 1

            #Make sure there are enough names for the vertices
            if vertices_count > len(vnames):
                raise ValueError('More vertices than vnames for '+topology_val)

            #Make a coordination number dictionary
            cn[vnames[vertices_count-1]] = int(vert_val[1])

        #Get edge centers
        elif 'edge_center' in line.lower():
            edge_center_val = line.lower().split('edge_center')[-1].strip()
            edge_center_val = [float(i) for i in edge_center_val.split()]

            #Add 0.0 dummy coordinate for 2D topologies
            if len(edge_center_val) == 2:
                if topology_val in dummya_list:
                    edge_center_val.insert(0,0.0)
                elif topology_val in dummyb_list:
                    edge_center_val.insert(1,0.0)
                else:
                    edge_center_val.append(0.0)

            #Make sure there are [x,y,z] coordinates
            if len(edge_center_val) != 3:
                bad = True
                continue

            edges_center.append(edge_center_val)

        #Get edge endpoint
        elif 'edge' in line.lower():
            edge_val = line.lower().split('edge')[-1].strip()
            edge_val = [float(i) for i in edge_val.split()]

            if len(edge_val) == 4:
                edge_head_val = edge_val[0:2]
                edge_tail_val = edge_val[2:]
                if topology_val in dummya_list:
                    edge_head_val.insert(0,0.0)
                    edge_tail_val.insert(0,0.0)
                elif topology_val in dummyb_list:
                    edge_head_val.insert(1,0.0)
                    edge_tail_val.insert(1,0.0)
                else:
                    edge_head_val.append(0.0)
                    edge_tail_val.append(0.0)
            elif len(edge_val) == 6:
                edge_head_val = edge_val[0:3]
                edge_tail_val = edge_val[3:]
            else:
                bad = True
                continue

            #Make sure there are [x,y,z] coordinates
            if len(edge_head_val) != 3:
                bad = True
                continue            

            edges_head.append(edge_head_val)
            edges_tail.append(edge_tail_val)

        #Store results for topology
        elif line.lower() == 'end':

            if not consider_2D and not three_dim:
                continue

            #Skip weirdly formatted cgd entries
            if bad or len(cn) != len(vertices) or len(edges_head) != len(edges_center):
                warnings.warn('Error: skipping '+topology_val+' because it is not formatted properly in .cgd file',Warning)
                continue            

            topologies_all.append(topology_val)
            groups_all.append(group_val)
            cellpars_all.append(cell_val)
            vertices_all.append(vertices)
            edges_center_all.append(edges_center)
            edges_head_all.append(edges_head)
            edges_tail_all.append(edges_tail)
            cn_all.append(cn)
            is_threedim_all.append(three_dim)

        #Ignore NC nets (assumed to be at bottom of .cgd file)
        elif 'nc nets' in line.lower():
            break

#Make folders to store topology CIFs
if not os.path.exists('template_database'):
    os.mkdir('template_database')
if not os.path.exists('template_errors'):
    os.mkdir('template_errors')
if not os.path.exists('template_2D_database') and consider_2D:
    os.mkdir('template_2D_database')

#Cycle through all topologies and make CIFs
for i in range(0,len(topologies_all)):

    #Flag for skipping CIF generation
    bad = False

    #Get all .cgd info for given topology, i
    topology = topologies_all[i]
    group = groups_all[i]
    cellpars = cellpars_all[i]
    vertices = vertices_all[i]
    edges_center = edges_center_all[i]
    edges_head = edges_head_all[i]
    edges_tail = edges_tail_all[i]
    cn_vec = cn_all[i]
    threedim = is_threedim_all[i]

    #Make list of vertex and edge center symbols
    sym_vertices = []
    for j in range(len(vertices)):
        sym_vertices.append(vnames[j])
    sym_collection = sym_vertices+[edge_center_name]*len(edges_center)

    #Get lattice vectors (using ASE function because it's easy)
    lattice_vectors = cellpar_to_cell(cellpars)

    #Get vertex and edge positions
    basis_collection = np.array(vertices+edges_center)

    #Make pymatgen structure
    try:
        pm_structure = pm.Structure.from_spacegroup(group,lattice_vectors,sym_collection,basis_collection)
    except:
        warnings.warn('Error: '+topology+'. Incompatible spacegroup and lattice constants',Warning)
        pm_structure.to(filename=os.path.join('template_errors',topology+'.cif'))
        continue        
    pm_structure.merge_sites(mode='delete')    

    #Calculate distance between edge centers and edge ends
    bd_list = []
    for j, edge_center_pos in enumerate(edges_center):
        n_edge_type = len(pm.Structure.from_spacegroup(group,lattice_vectors,[edge_center_name],[edge_center_pos]))
        dummy_edge = pm.Structure(lattice_vectors,[edge_center_name],[edges_head[j]])[0]
        dummy_center = pm.Structure(lattice_vectors,[edge_center_name],[edge_center_pos])[0]
        bd_list.extend([2*dummy_center.distance(dummy_edge)]*n_edge_type)

    _, unique_indices = np.unique(bd_list, return_index=True)
    unique_bond_dists = np.array(bd_list)[np.sort(unique_indices)].tolist()
    if np.abs(np.max(unique_bond_dists)-np.min(unique_bond_dists)) < tol:
        unique_bond_dists = [np.average(unique_bond_dists)]

    #Make lattice constants > bond dist
    if np.max(unique_bond_dists) < scale:
        extend = tol+scale
    else:
        extend = tol+np.max(unique_bond_dists)
    n_supercells = [np.ceil(extend/cellpars[0]),np.ceil(extend/cellpars[1]),np.ceil(extend/cellpars[2])]
    if n_supercells != [1,1,1]:
        pm_structure.make_supercell(n_supercells)

    #Get atoms of edge centers and vertices
    vertices_indices = [atom_idx for atom_idx, atom in enumerate(pm_structure) if atom.species_string != edge_center_name]
    edge_center_indices = [atom_idx for atom_idx, atom in enumerate(pm_structure) if atom.species_string == edge_center_name]

    #Make text for top of CIF
    top_text = 'data_'+topology+'\n'+'_audit_creation_date              '+datetime.today().strftime('%Y-%m-%d')+'\n'+"_audit_creation_method            'Pymatgen'\n"+"_symmetry_space_group_name_H-M    'P1'\n"+'_symmetry_Int_Tables_number       1\n'
    cellpar_text = 'loop_\n_symmetry_equiv_pos_as_xyz\n  x,y,z\n'+'_cell_length_a                    '+str(np.round(pm_structure.lattice.abc[0],4))+'\n'+'_cell_length_b                    '+str(np.round(pm_structure.lattice.abc[1],4))+'\n'+'_cell_length_c                    '+str(np.round(pm_structure.lattice.abc[2],4))+'\n'+'_cell_angle_alpha                 '+str(np.round(pm_structure.lattice.angles[0],4))+'\n'+'_cell_angle_beta                 '+str(np.round(pm_structure.lattice.angles[1],4))+'\n'+'_cell_angle_gamma                 '+str(np.round(pm_structure.lattice.angles[2],4))+'\n'
    pos_text = 'loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n'

    #Initialization
    bonded_pairs = [] #list for the indices of bonded vertices
    bonded_edge_centers = [] #list for the indices of bonded edge centers
    img_list = [] #list of image displacements
    d_list = [] #list of bond distances
    bonded_set_all = [] #list of all bonded sets
    bonded_edges_all = [] #list of all edges involved in bonds

    #Cycle through every vertex to find its bonded atoms
    for j, vertex_idx in enumerate(vertices_indices):

        #Initialization
        vertex_atom = pm_structure[vertex_idx] #Site object
        cn = cn_vec[vertex_atom.species_string] #int
        pm_structure[vertex_idx].index = j #store the index, excluding edge centers

        #Make string containing fract position for atom j
        pos_text += vertex_atom.species_string+str(j+1)+'     '+vertex_atom.species_string+'     '+str(np.round(vertex_atom.frac_coords[0],4))+'   '+str(np.round(vertex_atom.frac_coords[1],4))+'   '+str(np.round(vertex_atom.frac_coords[2],4))+'\n'

        #Find all edge centers connected to vertex j
        edge_overlap_indices = []
        for bond_dist in unique_bond_dists:
            edges_shell_temp = pm_structure.get_neighbors_in_shell(pm_structure[vertex_idx].coords,bond_dist/2,tol,include_index=True)
            edges_shell = [k for k in edges_shell_temp if k[0].species_string == edge_center_name]
            for edge_shell in edges_shell:
                if edge_shell[2] not in edge_overlap_indices:
                    edge_overlap_indices.append(edge_shell[2])
        edge_overlap_indices.sort()

        #Make sure the right number of edges are detected
        if len(edge_overlap_indices) != cn:
            warnings.warn('Error: '+topology+'. Incorrect number of edges',Warning)
            pm_structure.to(filename=os.path.join('template_errors',topology+'.cif'))
            bad = True
            break

        #Generate bond info
        vertex_overlap_indices = []
        bonded_set = []
        for bond_dist in unique_bond_dists:

            big_vertices_temp = pm_structure.get_neighbors_in_shell(pm_structure[vertex_idx].coords,bond_dist,2*tol,include_index=True,include_image=True)
            big_vertices = [k for k in big_vertices_temp if k[0].species_string != edge_center_name and k[2] != vertex_idx]
            big_vertices_indices = [b[2] for b in big_vertices]

            #Cycle through every edge center connected to vertex j
            for edge2_index in edge_overlap_indices:

                edge_overlap_atom = pm_structure[edge2_index]

                #Find bonded vertex connecting j and edge center
                vertices_shell_temp = pm_structure.get_neighbors_in_shell(edge_overlap_atom.coords,bond_dist/2,tol,include_index=True)    
                vertices_shell = [k for k in vertices_shell_temp if k[0].species_string != edge_center_name and k[2] in big_vertices_indices]
                if len(vertices_shell) != 1:
                    continue
                bonded_vertex_idx = vertices_shell[0][2]

                #Get the image of the bond
                possible_images = [k for k in big_vertices if k[2] == bonded_vertex_idx]
                for possible_image in possible_images:
                    dummy_pos = vertex_atom.coords+(possible_image[0].coords-vertex_atom.coords)/2
                    dummy_atom = pm.Structure(pm_structure.lattice.matrix,[edge_center_name],[dummy_pos],coords_are_cartesian=True)
                    if edge_overlap_atom.is_periodic_image(dummy_atom[0],tolerance=2*tol):
                        img_temp = possible_image[3].tolist()
                        d_img = possible_image[1]
                        z = [vertex_idx,bonded_vertex_idx,img_temp]
                        break
                if z in bonded_set_all:
                    continue

                #Store results
                bonded_set.append(z)
                bonded_set_all.append(z)
                vertex_overlap_indices.append(bonded_vertex_idx)
                img_list.append([int(ii) for ii in img_temp])
                d_list.append(d_img)

        #Check coordination number
        if len(vertex_overlap_indices) != cn:
            warnings.warn('Error: '+topology+'. Incorrect number of bonded vertices',Warning)
            pm_structure.to(filename=os.path.join('template_errors',topology+'.cif'))
            bad = True
            break

        #Add set of bonded pair of indices to list
        for vertex_overlap_idx in vertex_overlap_indices:
            bonded_pairs.append([vertex_idx,vertex_overlap_idx])
        
    if bad:
        continue

    #Make the bonding text for the CIF
    bond_text = 'loop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_distance\n_geom_bond_site_symmetry_2\n_ccdc_geom_bond_type\n'

    #For every bonded pair, get bonding/symmetry info
    done_dot_indices = [] #completed bond pairs with . symmetry
    bond_counts = [0]*len(vertices_indices)
    for j, bonded_pair in enumerate(bonded_pairs):

        #Get distance/image properties
        atom1 = pm_structure[bonded_pair[0]] #Vertex1
        atom2 = pm_structure[bonded_pair[1]] #Vertex2
        output_indices = [atom1.index+1,atom2.index+1] #indices to write in CIF
        img = img_list[j]
        d = d_list[j]

        #Make symmetry text
        if img == [0,0,0]:
            symmetry_sym = '.'
        else:
            symmetry_sym = '1_'+str(img[0]+5)+str(img[1]+5)+str(img[2]+5)

        #Complete bond text string
        if symmetry_sym == '.' and (output_indices in done_dot_indices or [output_indices[1],output_indices[0]] in done_dot_indices):
            continue
        bond_text += atom1.species_string+str(output_indices[0])+'     '+atom2.species_string+str(output_indices[1])+'    '+str(np.round(d,3))+'   '+symmetry_sym+'     S\n'
        if symmetry_sym == '.':
            done_dot_indices.append(output_indices)
            bond_counts[atom1.index] += 1
            bond_counts[atom2.index] += 1
        else:
            bond_counts[atom1.index] += 0.5
            bond_counts[atom2.index] += 0.5
            if [bonded_pair[1],bonded_pair[0],[-ii for ii in img]] not in bonded_set_all:
                warnings.warn('Error: '+topology+'. Missing symmetry counterpoint',Warning)
                pm_structure.to(filename=os.path.join('template_errors',topology+'.cif'))
                bad = True
                break                

    if bad:
        continue

    for j, bond_count in enumerate(bond_counts):
        if bond_count != cn_vec[pm_structure[vertices_indices[j]].species_string]:
            warnings.warn('Error: '+topology+'. Incorrect number of bonded vertices in symmetry flags',Warning)
            pm_structure.to(filename=os.path.join('template_errors',topology+'.cif'))
            bad = True
            break

    if bad:
        continue

    #Write the topology CIF
    if threedim:
        with open(os.path.join('template_database',topology+'.cif'),'w') as w:
            w.write(top_text+cellpar_text+pos_text+bond_text)
    else:
        with open(os.path.join('template_2D_database',topology+'.cif'),'w') as w:
            w.write(top_text+cellpar_text+pos_text+bond_text)        
    print('Success: '+topology)
    