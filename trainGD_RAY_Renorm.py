import scipy.io
from scipy.optimize import minimize
import numpy as np
from scipy.integrate import quad
#import numpy.random as np.random
np.random.seed(2018)

def trainGD_RAY(seq,eps):

	#eps = np.finfo(float).eps

	gamma_0 = np.random.rand()
	eta_0 = np.random.rand() #2*gamma_0
	mu_0 = np.random.rand()

	# input_data = scipy.io.loadmat('4Kern_newSNS_10seqT100000_highMuandfreq0.15.mat')
	# seq = input_data['Seq2'][1][0][0]

	# seq = seq[:300]

	T = seq[-1]#-seq[0]
	Delta = len(seq)/T

	#seq = np.array([1.,2.,3.,4.,5.])

	bnds = ((0,None),(0,None),(0,None))

	def logGD_RAY(RAY_coeffs):

		def funcRAY(x,gamma,eta):
			return gamma*x*np.exp(-eta*np.power(x,2))

		gamma = RAY_coeffs[1];

		eta = RAY_coeffs[2];

		mu = RAY_coeffs[0]

		if (gamma/(2*eta) < 1.) and (gamma/(2*eta) >= 0.):
			mu = mu#(1.-gamma/(2*eta))*Delta;
		else:
			mu = mu#0. 
			#return np.inf

		intens = np.zeros(len(seq));

		compens = mu*T;

		for i in range(0,len(seq)):

			intens[i] += mu;

			#############  REVER ESSA PARTE !!!!!!!!!!!!! ##################################################

			compens += (gamma/(2*eta))*(1-np.exp(-eta*np.power(T-seq[i],2)))#quad(funcRAY,0,T-seq[i], args=(gamma,eta))[0]

			print('compens: '+repr(compens))

			for j in range(0,i):

				intens[i] += gamma*(seq[i]-seq[j])*np.exp(-eta*np.power(seq[i] - seq[j],2))			

			#print('intens_i: '+repr(intens))

		intens[intens < 0.] = 0.	

		print ('Loglikelihood Train GD: ' + repr(np.sum(np.nan_to_num(np.log(intens))) - compens) + '\n')

		return - np.sum(np.nan_to_num(np.log(intens))) + max(compens,0.)

	par = minimize(logGD_RAY, [mu_0, gamma_0, eta_0], method='Nelder-Mead', tol=1e-2, options={'maxiter':10})

	print('Final Parameters: '+ repr(par.x)+'\n')

	RAY_statcriter = par.x[1]/(2*par.x[2])

	print('RAY_statcriter: ' + repr(RAY_statcriter))

	fin_llh = logGD_RAY(par.x)

	fin_llh = (-1)*fin_llh

#	if np.isinf(fin_llh) or np.isnan(fin_llh) or (fin_llh > 0.):

	par_renorm_gamma = [Delta*(1.-1./(1.+eps)),par.x[1]/(RAY_statcriter*(1+eps)),par.x[2]]

	par_renorm_eta = [Delta*(1.-1./(1.+eps)),par.x[1],par.x[2]*(RAY_statcriter*(1+eps))]

	par_renorm_sqrt = [Delta*(1.-1./(1.+eps)),par.x[1]/np.sqrt(RAY_statcriter*(1+eps)),par.x[2]*np.sqrt(RAY_statcriter*(1+eps))]

	llh_renorm_gamma = logGD_RAY(par_renorm_gamma)

	llh_renorm_eta = logGD_RAY(par_renorm_eta)

	llh_renorm_sqrt = logGD_RAY(par_renorm_sqrt)

	llh_renorm_gamma *= -1

	llh_renorm_eta *= -1

	llh_renorm_sqrt  *= -1

	print('par_renorm_gamma: '+repr(par_renorm_gamma))

	print('par_renorm_eta: '+repr(par_renorm_eta))

	print('par_renorm_sqrt: '+repr(par_renorm_sqrt))

	print('llh_renorm_gamma: '+repr(llh_renorm_gamma))

	print('llh_renorm_eta: '+repr(llh_renorm_eta))

	print('llh_renorm_sqrt: '+repr(llh_renorm_sqrt))

	K1_Param = {'RAY_coeffs': par.x, 'K1_Type': 'RAY', 'RAY_statcriter': par.x[1]/par.x[2], 'final_llh': fin_llh, 'par_renorm_gamma': par_renorm_gamma,\
	'llh_renorm_gamma': llh_renorm_gamma, 'par_renorm_eta': par_renorm_eta, 'llh_renorm_eta': llh_renorm_eta, 'par_renorm_sqrt': par_renorm_sqrt, 'llh_renorm_sqrt':llh_renorm_sqrt}

	# else:

	# 	llh_renorm_gamma = fin_llh

	# 	llh_renorm_eta = fin_llh

	# 	llh_renorm_sqrt = fin_llh

	# 	K1_Param = {'RAY_coeffs': par.x, 'K1_Type': 'RAY', 'RAY_statcriter': par.x[1]/par.x[2], 'final_llh': fin_llh,\
	# 	'llh_renorm_gamma': llh_renorm_gamma, 'llh_renorm_eta': llh_renorm_eta, 'llh_renorm_sqrt':llh_renorm_sqrt}

	return K1_Param
