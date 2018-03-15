import scipy.io
from scipy.optimize import minimize
import numpy as np
from scipy.integrate import quad
#import numpy.random as np.random
np.random.seed(2018)

def trainGD_EXP(seq,eps):

	#eps = np.finfo(float).eps

	alpha_0 = np.random.rand()
	beta_0 = np.random.rand() #2*alpha_0
	mu_0 = np.random.rand()

	# input_data = scipy.io.loadmat('4Kern_newSNS_10seqT100000_highMuandfreq0.15.mat')
	# seq = input_data['Seq2'][1][0][0]

	# seq = seq[:300]

	T = seq[-1]#-seq[0]
	Delta = len(seq)/T

	#seq = np.array([1.,2.,3.,4.,5.])

	bnds = ((0,None),(0,None),(0,None))

	def logGD_EXP(EXP_coeffs):

		def funcexp(x,alpha,beta):
			return alpha*np.exp(-beta*x)

		alpha = EXP_coeffs[1];

		beta = EXP_coeffs[2];

		mu = EXP_coeffs[0]

		if (alpha/beta < 1.) and (alpha/beta >= 0.):
			mu = mu#(1.-alpha/beta)*Delta;
		else:
			mu = mu#0. 
			#return np.inf

		intens = np.zeros(len(seq));

		compens = mu*T;

		for i in range(0,len(seq)):

			intens[i] += mu;

			compens += (alpha/beta)*(1-np.exp(-beta*(T-seq[i])))#quad(funcexp,0,T-seq[i], args=(alpha,beta))[0]

			for j in range(0,i):

				intens[i] += alpha*np.exp(-beta*(seq[i] - seq[j]))			

		print ('Loglikelihood Train GD: ' + repr(np.sum(np.nan_to_num(np.log(intens))) - compens) + '\n')

		return - np.sum(np.nan_to_num(np.log(intens))) + compens

	par = minimize(logGD_EXP, [mu_0, alpha_0, beta_0], method='Nelder-Mead', tol=1e-2, options={'maxiter':10})

	print('Final Parameters: '+ repr(par.x)+'\n')

	EXP_statcriter = par.x[1]/par.x[2]

	print('EXP_statcriter:' + repr(EXP_statcriter))

	fin_llh = logGD_EXP(par.x)

	fin_llh = (-1)*fin_llh

#	if np.isinf(fin_llh) or np.isnan(fin_llh) or (fin_llh > 0.):

	par_renorm_alpha = [Delta*(1.-1./(1.+eps)),par.x[1]/(EXP_statcriter*(1+eps)),par.x[2]]

	par_renorm_beta = [Delta*(1.-1./(1.+eps)),par.x[1],par.x[2]*(1+eps)]

	par_renorm_sqrt = [Delta*(1.-1./(1.+eps)),par.x[1]/np.sqrt(EXP_statcriter*(1+eps)),par.x[2]*np.sqrt(1+eps)]

	print('par_renorm_alpha: '+repr(par_renorm_alpha))

	print('par_renorm_beta: ' + repr(par_renorm_beta))

	print('par_renorm_sqrt: ' + repr(par_renorm_sqrt))

	llh_renorm_alpha = logGD_EXP(par_renorm_alpha)

	llh_renorm_beta = logGD_EXP(par_renorm_beta)

	llh_renorm_sqrt = logGD_EXP(par_renorm_sqrt)

	llh_renorm_alpha *= -1

	llh_renorm_beta *= -1

	llh_renorm_sqrt *= -1

	print('llh_renorm_alpha: '+repr(llh_renorm_alpha))

	print('llh_renorm_beta: '+repr(llh_renorm_beta))

	print('llh_renorm_sqrt: '+ repr(llh_renorm_sqrt))


#	else:

		# llh_renorm_alpha = fin_llh

		# llh_renorm_beta = fin_llh

		# llh_renorm_sqrt = fin_llh

	K1_Param = {'EXP_coeffs': par.x, 'K1_Type': 'EXP', 'EXP_statcriter': par.x[1]/par.x[2], 'final_llh': fin_llh, 'llh_renorm_alpha': llh_renorm_alpha, 		'llh_renorm_beta': llh_renorm_beta, 'llh_renorm_sqrt':llh_renorm_sqrt}

	return K1_Param

