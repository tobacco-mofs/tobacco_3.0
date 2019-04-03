import numpy as np
from scipy.optimize import minimize, basinhopping

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

def scale(all_SBU_coords,a,b,c,ang_alpha,ang_beta,ang_gamma,max_le,num_vertices,Bstar,alpha,num_edges,PATIENCE):

	max_length = 0
	for line in all_SBU_coords:
		for length in [np.linalg.norm(s[1]) for s in line[1]]:
			if length > max_length:
				max_length = length

	scale_guess = (max_length / max_le)
	#scale_guess = 1.0
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
	x_bounds = tuple([(None,None) for x in covars_values])

	Bstar_inv = np.linalg.inv(Bstar)

	print 'scaling unit cell and vertex positions...'
	print ''

	if PATIENCE:
		res = basinhopping(objective, init_variables, niter=100, T=1.0, stepsize=0.5, 
				 	   minimizer_kwargs={'args':(ncra,ncca,alpha,num_edges,num_vertices,Bstar_inv,all_SBU_ip), 
				 				   'method':'L-BFGS-B', 
				 				   'bounds':((0,None),(0,None),(0,None),(20,160),(20,160),(20,160)) + x_bounds})
	else:
		res = minimize(objective, init_variables, args=(ncra,ncca,alpha,num_edges,num_vertices,Bstar_inv,all_SBU_ip),
						method='L-BFGS-B',
						bounds=((0,None),(0,None),(0,None),(20,160),(20,160),(20,160)) + x_bounds,
						options={'disp':False, 'gtol':1E-12, 'ftol':1E-12})

	sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = res.x[0:6]
	sc_covar = res.x[6:].reshape(ncra,ncca)

	return(sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma,sc_covar,Bstar_inv,max_length)

