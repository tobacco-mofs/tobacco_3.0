import numpy as np
import networkx as nx
import math
import re
from scaled_embedding2coords import omega2coords
from write_cifs import PBC3DF, PBC3DF_sym

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import animation
from matplotlib.lines import Line2D

def nn(string):
	return re.sub('[^a-zA-Z]','', string)

def nl(string):
	return re.sub('[^0-9]','', string)

def roundup(x):
	return int(math.ceil(x / 10.0)) * 10

pi = np.pi

def scaling_callback_animation(callbackresults, alpha, Bstar_inv, ncra, ncca, num_vertices, num_edges, TG, template, g, cfiles):
					
	frames = []
	frames_append = frames.append

	for l in callbackresults:
		
		X, fX = l
		uc_params = X[0:6]
		rc_a, rc_b, rc_c, rc_alpha, rc_beta, rc_gamma = uc_params
		rc_covars = np.array(X[6:])
		rc_covars = rc_covars.reshape(ncra,ncca)
		rc_Alpha = np.r_[alpha[0:num_edges-num_vertices+1,:], rc_covars]
		rc_omega_plus = np.dot(Bstar_inv, rc_Alpha)
		rc_coords = omega2coords(TG, rc_omega_plus, uc_params, num_vertices, template, g, cfiles)

		ax_rc = rc_a
		ay_rc = 0.0
		az_rc = 0.0
		bx_rc = rc_b * np.cos(rc_gamma * pi/180.0)
		by_rc = rc_b * np.sin(rc_gamma * pi/180.0)
		bz_rc = 0.0
		cx_rc = rc_c * np.cos(rc_beta * pi/180.0)
		cy_rc = (rc_c * rc_b * np.cos(rc_alpha * pi/180.0) - bx_rc * cx_rc) / by_rc
		cz_rc = (rc_c ** 2.0 - cx_rc ** 2.0 - cy_rc ** 2.0) ** 0.5
		rc_unit_cell = np.asarray([[ax_rc,ay_rc,az_rc],[bx_rc,by_rc,bz_rc],[cx_rc,cy_rc,cz_rc]]).T

		frames_append([rc_coords, rc_unit_cell])

	return frames

def write_scaling_callback_animation(frames, prefix):

	G = nx.Graph()
	coord_frames = []
	coord_frames_append = coord_frames.append
	inds = []
	inds_append = inds.append
	norm_coord = np.array([0,0,0])
	norm_coord_dict = {}
	fcoord_dict = {}

	for l in frames[0][0]:
		name, cif, fcoords, edges = l
		inds_append(int(nl(name)))

	max_ind = max(inds)

	for l in frames[0][0]:

		name, cif, fcoords, edges = l
		norm_coord_dict[name] = fcoords
		fcoord_dict[name] = np.array(fcoords)

	transform_dict = dict((k,[]) for k in fcoord_dict)

	for l in frames[0][0]:

		name, cif, fcoords, edges = l

		for e in edges:
			a1,a2 = e[1]
			dist, sym = PBC3DF_sym(fcoord_dict[a1], fcoord_dict[a2])
			if sym == '.':
				G.add_edge(*e[1])
				transform_dict[a1].append((a1, (0,0,0)))
				transform_dict[a2].append((a2, (0,0,0)))
			else:
				sym = sym.split('_')[1]
				sym = np.array([float(sym[0]),float(sym[1]),float(sym[2])])
				sym = sym - 5.0
				max_ind += 1
				na = nn(a2) + str(max_ind)
				G.add_edge(a1, na)
				transform_dict[a2].append((na, (sym[0],sym[1],sym[2])))

		duplicate_vecs = []
		for dim in (0,1,2):
			init = [0,0,0] 
			if abs(fcoords[dim] - 1.0) < 0.05:
				init[dim] -= 1.0
			elif abs(fcoords[dim]    ) < 0.05:
				init[dim] += 1.0
			duplicate_vecs.append(init)

		for dv in duplicate_vecs:
			dup_vec = np.array(dv)
			if np.any(dup_vec):
				max_ind += 1
				new_name = nn(name) + str(max_ind)
				transform_dict[name].append((new_name, (dv[0], dv[1], dv[2])))

				for e in edges:
					a1,a2 = e[1]
					dist, sym = PBC3DF_sym(fcoords + dup_vec, fcoord_dict[a2])
					if sym == '.':
						G.add_edge(new_name, a2)

	for a in transform_dict:
		transform_dict[a] = set(transform_dict[a])

	for f in frames:

		coords = []
		coords_append = coords.append

		fc, rc_unit_cell = f
		fcoord_dict = dict((s[0],np.array(s[2])) for s in fc)

		for l in fc:

			name, cif, fcoords, edges = l
			fcoords = PBC3DF(norm_coord_dict[name], fcoords)
			ccoords = np.dot(rc_unit_cell, fcoords) 
			coords_append([name, ccoords ])

			if len(transform_dict[name]) > 1:
				for na in transform_dict[name]:
					new_name, tvec = na
					if new_name != name:
						tvec = np.asarray(tvec)
						new_ccoords = np.dot(rc_unit_cell, fcoords + tvec) 
						coords_append([new_name, new_ccoords])

		coord_frames_append(coords)

	for i in range(len(coord_frames)):
		f = coord_frames[i]
		coord_frames[i] = sorted(f, key = lambda x : int(nl(x[0])))

	with open(prefix + '_scale_animation.xyz', 'w') as out:
		fc = 0
		for f in coord_frames:
			fc += 1
			out.write(str(len(f)) + '\n')
			out.write('frame: ' + str(fc) + '\n')
			for l in f:
				name, coord = l
				x,y,z = coord
				out.write('{:5} {:>10.5f} {:>10.5f} {:>10.5f}'.format(name, x, y, z))
				out.write('\n')

	with open(prefix + '_add_bonds.tcl', 'w') as out:
			for b in G.edges():
				ind0 = int(nl(b[0])) - 1
				ind1 = int(nl(b[1])) - 1
				out.write('topo addbond ' + str(ind0) + ' ' + str(ind1) + '\n')

def animate_objective_minimization(callbackresults, prefix, bitrate=1800, fps=30, font_size=16, time=5.0):

	plt.rcParams.update({'font.size': font_size})

	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=bitrate)

	uc_params = []
	uc_params_append = uc_params.append
	fX_vals = []
	fX_vals_append = fX_vals.append
	res_vals = []
	res_vals_append = res_vals.append

	max_len = 0
	max_angle = 0

	for l in callbackresults:

		X, fX = l
		a,b,c,alpha,beta,gamma = X[0:6]

		max_len_current = max((a,b,c))
		if max_len_current > max_len:
			max_len = max_len_current
		max_angle_current = max((alpha,beta,gamma))
		if max_angle_current > max_angle:
			max_angle = max_angle_current

		res = X[6:]
		avg_res = np.average(res)
		uc_params.append([a, b, c, alpha, beta, gamma])
		fX_vals.append(fX)
		res_vals.append(avg_res)

	iterations = np.array(range(len(fX_vals)))
	uc_params = np.asarray(uc_params)
	fX_vals = np.asarray(fX_vals)/max(fX_vals)
	res_vals = np.asarray(res_vals)

	nframes = len(iterations)
	time *= 1000
	delay = time/nframes

	class SubplotAnimation(animation.TimedAnimation):
		def __init__(self):
			
			fig = plt.figure(figsize=(16.0,5.5))
			matplotlib.pyplot.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.25, hspace=None)

			ax1 = fig.add_subplot(1, 3, 1)
			ax2 = fig.add_subplot(1, 3, 2)
			ax3 = fig.add_subplot(1, 3, 3)

			self.t     = iterations
			self.x     = iterations
			self.fX    = fX_vals
			self.a     = uc_params[:,0]
			self.b     = uc_params[:,1]
			self.c     = uc_params[:,2]
			self.alpha = uc_params[:,3]
			self.beta  = uc_params[:,4]
			self.gamma = uc_params[:,5]

			ax1.set_xlabel('Iterations')
			ax1.set_ylabel('Normalized Objective Function')
			self.line11 = Line2D([], [], color='black', linewidth=2)
			ax1.add_line(self.line11)
			ax1.set_xlim(0, max(iterations))
			ax1.set_ylim(0, 0.10)
			ax1.xaxis.set_major_formatter(plt.NullFormatter())

			ax2.set_xlabel('Iterations')
			ax2.set_ylabel('Unit Cell Lengths / ' + r'$\AA$')
			self.line21 = Line2D([], [], color='black', linewidth=2)
			self.line22 = Line2D([], [], color='blue' , linewidth=2)
			self.line23 = Line2D([], [], color='red'  , linewidth=2)
			ax2.add_line(self.line21)
			ax2.add_line(self.line22)
			ax2.add_line(self.line23)
			ax2.set_xlim(0, max(iterations))
			ax2.set_ylim(5, max_len + 5.0)
			ax2.xaxis.set_major_formatter(plt.NullFormatter())

			ax3.set_xlabel('Iterations')
			ax3.set_ylabel('Unit Cell Angles / degrees')
			self.line31 = Line2D([], [], color='black', linewidth=2)
			self.line32 = Line2D([], [], color='blue' , linewidth=2)
			self.line33 = Line2D([], [], color='red'  , linewidth=2)
			ax3.add_line(self.line31)
			ax3.add_line(self.line32)
			ax3.add_line(self.line33)
			ax3.set_xlim(0 , max(iterations))
			ax3.set_ylim(60, max_angle + 10.0)
			ax3.xaxis.set_major_formatter(plt.NullFormatter())

			animation.TimedAnimation.__init__(self, fig, interval=delay, blit=True)

		def _draw_frame(self, framedata):
			
			i = framedata

			self.line11.set_data(self.x[:i], self.fX[:i])

			self.line21.set_data(self.x[:i], self.a[:i])
			self.line22.set_data(self.x[:i], self.b[:i])
			self.line23.set_data(self.x[:i], self.c[:i])

			self.line31.set_data(self.x[:i], self.alpha[:i])
			self.line32.set_data(self.x[:i], self.beta[:i])
			self.line33.set_data(self.x[:i], self.gamma[:i])

			self._drawn_artists = [self.line11, self.line21, self.line22,
								   self.line23, self.line31, self.line32,
								   self.line33]

		def new_frame_seq(self):
			return iter(range(self.t.size))

		def _init_draw(self):
			lines = [self.line11, self.line21, self.line22,
					 self.line23, self.line31, self.line32,
					 self.line33]
			for l in lines:
				l.set_data([], [])

	ani = SubplotAnimation()
	ani.save(prefix + '_ts.mp4')

