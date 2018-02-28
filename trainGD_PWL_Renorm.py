import scipy.io
from scipy.optimize import minimize
import numpy as np
from scipy.integrate import quad
from scipy.special import lambertw
#import numpy.random as np.random

def trainGD_PWL(seq,eps):

	T = seq[-1]-seq[0]
	Delta = len(seq)/T

	mu_0 = Delta
	# K_0 = 0.59571517
	# c_0 = 1.94173945
	# p_0 = 2.13085285

	K_0 = np.random.rand()
	c_0 = 2*K_0
	p_0 = 1. + np.random.rand()

	# input_data = scipy.io.loadmat('4Kern_newSNS_10seqT100000_highMuandfreq0.15.mat')
	# seq = input_data['Seq2'][1][0][0]

	# seq = seq[:300]

	#seq = np.array([1.,2.,3.,4.,5.])

	bnds = ((0,None),(0,None),(0,None))

	def logGD_PWL(PWL_coeffs):

		def funcpwl(x,K,c,p):
			return K*np.power(x+c,-p)

		K = PWL_coeffs[0];

		c = PWL_coeffs[1];

		p = PWL_coeffs[2];

		mu = PWL_coeffs[3]
		
		phi = K*np.power(c,1.-p)/(1.-p)

		if (phi < 1.) and (K >= 0.) and (c > 0.) and (p > 1.):
		 	mu = (1.-phi)*Delta;
		else:
			mu = 0. 
			return np.inf

		intens = np.zeros(len(seq));

		compens = mu*T;

		for i in range(0,len(seq)):

			intens[i] += mu;

			compens += K*np.power(T-seq[i]+c,1-p)/(1-p) - K*np.power(c,1-p)/(1-p);#quad(funcpwl,0,T-seq[i], args=(K,c,p))[0] #(alpha/beta)*(1-np.exp(-beta*(T-seq[i])))

			for j in range(0,i):

				intens[i] += K*np.power((seq[i] - seq[j])+c,-p)			

		print ('Loglikelihood Train GD: ' + repr(np.sum(np.nan_to_num(np.log(intens))) - compens) + '\n')

		return - np.sum(np.nan_to_num(np.log(intens))) + compens

	par = minimize(logGD_PWL, [K_0, c_0, p_0, mu_0], method='nelder-mead', tol=1e-2, options={'maxiter':10})

	print('Final Parameters: '+ repr(par.x)+'\n')

	PWL_statcriter = ((par.x[2]-1)*par.x[0]*(par.x[1]**(1-par.x[2])))

	fin_llh = logGD_PWL(par.x)

	fin_llh = (-1)*fin_llh

	if np.isinf(fin_llh):

		par_renorm_K = [par.x[0]/(PWL_statcriter*(1+eps)),par.x[1],par.x[2],par.x[3]]

		par_renorm_c = [par.x[0],np.power((1+eps)*PWL_statcriter,1/par.x[2])*par.x[1],par.x[2],par.x[3]]

		Delta_p = PWL_statcriter*(1+eps)*(par.x[2]-1)*np.power(par.x[1],par.x[2]-1.)

		par_renorm_p = [par.x[0],par.x[1],1+ lambertw(Delta_p*np.log(par.x[1])).real/np.log(par.x[1]),par.x[3]]

		par_renorm_Kc = [par.x[0]/np.sqrt(PWL_statcriter*(1+eps)),par.x[1]*np.power(PWL_statcriter*(1+eps),1/(2*(par.x[2]-1))),par.x[2],par.x[3]]

		Delta_Kp = np.sqrt(PWL_statcriter*(1+eps))*(par.x[2]-1)*np.power(par.x[1],par.x[2]-1.)

		par_renorm_Kp = [par.x[0]/np.sqrt(PWL_statcriter*(1+eps)),par.x[1],1+lambertw(Delta_Kp*np.log(par.x[1])).real/np.log(par.x[1]),par.x[3]]

		# par_renorm_Kc = [par.x[0]/np.sqrt(PWL_statcriter*(1+eps)),np.power((1+eps)*PWL_statcriter,\
		# 	1/(2*par.x[2]))*par.x[1],par.x[2],par.x[3]]

		# par_renorm_Kp = [par.x[0]/np.sqrt(PWL_statcriter*(1+eps)),par.x[1],par.x[2]+ \
		# np.log(np.sqrt((1+eps)*PWL_statcriter))/np.log(par.x[1]),par.x[3]]

		# par_renorm_sqrt = [par.x[0]/(PWL_statcriter*np.power(1+eps,1/3)),np.power(np.power(1+eps,1/3)*PWL_statcriter,\
		# 	1/(par.x[2]+ np.log(np.power(1+eps,1/3)*PWL_statcriter)/np.log(np.power(np.power(1+eps,1/3)*PWL_statcriter,1/par.x[2])*par.x[1])))*par.x[1],\
		# par.x[2]+ np.log(np.power(1+eps,1/3)*PWL_statcriter)/np.log(np.power(np.power(1+eps,1/3)*PWL_statcriter,1/par.x[2])*par.x[1]),par.x[3]]

		# par_renorm_sqrt = [par.x[0]/(PWL_statcriter*np.power(1+eps,1/3)),par.x[1]*np.power(PWL_statcriter,1/par.x[2]+np.log(PWL_statcriter*\
		# 	np.power(1+eps,1/3))/np.log(par.x[1]))*np.power(1+eps,1/3*(par.x[2]+np.log(PWL_statcriter*np.power(1+eps,1/3))/np.log(par.x[1]))),\
		# par.x[2]+np.log(PWL_statcriter*np.power(1+eps,1/3))/np.log(par.x[1]),par.x[3]]

		llh_renorm_K = logGD_PWL(par_renorm_K)

		llh_renorm_c = logGD_PWL(par_renorm_c)

		llh_renorm_p = logGD_PWL(par_renorm_p)

		llh_renorm_Kc = logGD_PWL(par_renorm_Kc)

		llh_renorm_Kp = logGD_PWL(par_renorm_Kp)

		llh_renorm_K *= -1

		llh_renorm_c *= -1

		llh_renorm_p *= -1

		llh_renorm_Kc *= -1

		llh_renorm_Kp *= -1

	else:

		llh_renorm_K = fin_llh

		llh_renorm_c = fin_llh

		llh_renorm_p = fin_llh

		llh_renorm_Kc = fin_llh

		llh_renorm_Kp = fin_llh

	K1_Param = {'PWL_coeffs': par.x, 'K1_Type': 'PWL', 'PWL_statcriter': ((par.x[2]-1)*par.x[0]*(par.x[1]**(1-par.x[2]))), 'final_llh': fin_llh,\
	 'llh_renorm_K': llh_renorm_K, 'llh_renorm_c': llh_renorm_c, 'llh_renorm_p': llh_renorm_p, 'llh_renorm_Kc': llh_renorm_Kc,\
	  'llh_renorm_Kp': llh_renorm_Kp}#<1)*(par.x[2]>1)}

	return K1_Param
