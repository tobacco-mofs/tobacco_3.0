import numpy as np
from scipy.optimize import minimize, basinhopping
import random

pi = np.pi

def metric_tensor(l):
	a,b,c,alpha,beta,gamma = l
	Z = np.array([[a*a, a*b*np.cos(gamma * (pi/180)), a*c*np.cos(beta * (pi/180))],
			      [a*b*np.cos(gamma * (pi/180)), b*b, b*c*np.cos(alpha * (pi/180))],
			  	  [a*c*np.cos(beta * (pi/180)), b*c*np.cos(alpha * (pi/180)), c*c]])
	return Z

def objective(V, ncra, ncca, Alpha, ne, nv, Bstar_inv, SBU_IP):
	P = V[0:6]
	Z = metric_tensor(P)
	X = np.array(V[6:])
	X = X.reshape(ncra,ncca)
	var_alpha = np.r_[Alpha[0:ne-nv+1,:], X]

	omp = np.dot(Bstar_inv, var_alpha)
	g = np.dot(omp, np.dot(Z,omp.T))

	O1 = 0
	for sbu in SBU_IP:
		O2 = 0
		for i in sbu[1]:
			O2 = O2 + (i[2] - g[i[0] - 1][i[1] - 1])**2
		O1 = O1 + O2

	return O1

def scale(all_SBU_coords,a,b,c,ang_alpha,ang_beta,ang_gamma,max_le,num_vertices,Bstar,alpha,num_edges,FIX_UC,SCALING_ITERATIONS,PRE_SCALE,SCALING_CONVERGENCE_TOLERANCE,SCALING_STEP_SIZE):

	scale_tol = SCALING_CONVERGENCE_TOLERANCE
	scale_eps = SCALING_STEP_SIZE
	max_length = 0
	for line in all_SBU_coords:
		for length in [np.linalg.norm(s[1]) for s in line[1]]:
			if length > max_length:
				max_length = length

	if PRE_SCALE == 'none':
		scale_guess = 1.0
	else:
		scale_guess = (max_length / max_le) * PRE_SCALE

	all_SBU_ip = []
	all_SBU_ip_append = all_SBU_ip.append
	for sbu in all_SBU_coords:
		SBU_ip = []
		SBU_ip_append = SBU_ip.append
		w = len(sbu[1])
		for i in range(w):
			ivec = sbu[1][i][1]
			iind = sbu[1][i][0]
			for j in range(i, w):
				jvec = sbu[1][j][1]
				jind = sbu[1][j][0]
				dot  = np.dot(ivec,jvec)
				SBU_ip_append([iind,jind,dot])
		all_SBU_ip_append((sbu[0], SBU_ip))
	ncra = num_vertices - 1
	ncca = 3

	covars_values = []
	covars_values_append = covars_values.append
	for i in range(ncra):
		for j in range(ncca):
			covars_values_append(0)
	
	init_variables = [scale_guess * a,scale_guess * b, scale_guess * c,ang_alpha,ang_beta,ang_gamma] + covars_values

	if np.any(FIX_UC):

		uc_bounds = []
		uc_bounds_append = uc_bounds.append
		for f,p in zip(FIX_UC, init_variables[0:6]):
			if f:
				uc_bounds_append((p,p))
			else:
				uc_bounds_append((0,None))
		uc_bounds = tuple(uc_bounds)
	else:
		uc_bounds = ((0,None),(0,None),(0,None),(20,160),(20,160),(20,160))

	x_bounds = tuple([(None,None) for x in covars_values])
	bounds = uc_bounds + x_bounds

	Bstar_inv = np.linalg.inv(Bstar)

	print 'scaling unit cell and vertex positions...'
	print '' 

	niter = SCALING_ITERATIONS
	uc_press = 0.0001
	covars_perturb = 0.025
	callbackresults = [[init_variables, objective(init_variables,ncra,ncca,alpha,num_edges,num_vertices,Bstar_inv,all_SBU_ip)]]

	def callbackF(X):
		var_string = ''
		for val in X:
			var_string += str(val) + '   '
		callbackresults.append([X, objective(X,ncra,ncca,alpha,num_edges,num_vertices,Bstar_inv,all_SBU_ip)])

	for it in range(niter):

		res = minimize(objective, init_variables, args=(ncra,ncca,alpha,num_edges,num_vertices,Bstar_inv,all_SBU_ip),
						method='L-BFGS-B',
						bounds=bounds,
						options={'disp':False, 'gtol':scale_tol, 'ftol':scale_tol, 'eps':scale_eps},
						callback=callbackF)

		uc_params = res.x[0:6]
		uc_params = [i - (i * uc_press) for i in uc_params]
		mult = [random.choice([-1, 1]) for i in range(len(res.x[6:]))]
		covars = [i + j * (i * covars_perturb) for i,j in zip(res.x[6:], mult)]
		init_variables = uc_params + covars

	if niter != 0:	
		sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = res.x[0:6]
		sc_covar = res.x[6:].reshape(ncra,ncca)
	else:
		init_variables = np.asarray(init_variables)
		sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = init_variables[0:6]
		sc_covar = init_variables[6:].reshape(ncra,ncca)

	ab = [a/b, sc_a/sc_b]
	ac = [a/c, sc_a/sc_c]
	bc = [b/c, sc_b/sc_c]
	alpha = [ang_alpha, sc_alpha]
	beta = [ang_beta, sc_beta]
	gamma = [ang_gamma, sc_gamma]
	covar = [np.average(abs(np.array(callbackresults[0][0][6:]) - np.array(callbackresults[-1][0][6:])))]
	final_obj = [callbackresults[-1][1]]
	scaling_data = [ab, ac, bc, alpha, beta, gamma, covar, final_obj]

	return(sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma,sc_covar,Bstar_inv,max_length,callbackresults,ncra,ncca,scaling_data)

