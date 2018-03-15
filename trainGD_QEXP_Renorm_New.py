import scipy.io
from scipy.optimize import minimize
import numpy as np
from scipy.integrate import quad
import math
#import numpy.random as np.random
np.random.seed(2018)

def trainGD_QEXP(seq,eps):

	alpha_0 = np.random.rand()
	beta_0 = 2*alpha_0
	q_0 = 1. + np.random.rand()
	mu_0 = np.random.rand()

	epsilon = np.finfo(float).eps

	# input_data = scipy.io.loadmat('4Kern_newSNS_10seqT100000_highMuandfreq0.15.mat')
	# seq = input_data['Seq2'][1][0][0]

	# seq = seq[:300]

	T = seq[-1]#-seq[0]
	Delta = len(seq)/T

	#seq = np.array([1.,2.,3.,4.,5.])

	bnds = ((0,None),(0,None),(0,None))

	def logGD_QEXP(QEXP_coeffs):

		def funcqexp(x,alpha,q):

			print('q: '+ repr(q))

			if abs(q-1.) <= 0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

				return alpha*np.exp((-1)*x)

			elif (q != 1.):

				return np.max(0.,alpha*np.power(1+(q-1)*x+epsilon,1/(1-q)))

		alpha = QEXP_coeffs[1];

		q = QEXP_coeffs[2]

		mu = QEXP_coeffs[0]

		if abs(q-1.) <= 0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

			if (alpha < 1.) and (alpha >= 0.):

				mu = mu#(1.-alpha)*Delta;

			else:

				mu = mu#0. 
				#return np.inf

		#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
		elif (q != 1.):

			if (alpha/(2-q) < 1.) and (alpha/(2-q) > 0.) and (alpha > 0.):

				mu = mu#(1.- alpha/(2-q))*Delta;

			else:

				mu = mu#0.
				#return np.inf

		else:

			mu = Delta

		intens = np.zeros(len(seq));

		compens = mu*T;

		for i in range(0,len(seq)):

			intens[i] += mu;

			if abs(q-1.) <= 0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

				compens += alpha*(1-np.exp((-1)*(T-seq[i])))

			#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
			elif (q > 1.) and (q < 2.):

				compens += (alpha/(2-q))*(np.power(1+(q-1)*(T-seq[i])+epsilon, (2-q)/(1-q))-1)

			elif (q < 1.) and (T-seq[i]>=1/(1-q)):

				compens += alpha/(2-q)

			elif (q < 1.) and (T - seq[i] < 1/(1-q)):

				compens += (alpha/(q-2.))*np.power(1.+(q-1.)*(T-seq[i])+epsilon,(2-q)/(1-q)) - alpha/(q-2)

			else:  

			#compens += (alpha/beta)*(1-np.exp(-beta*(T-seq[i])))#quad(funcexp,0,T-seq[i], args=(alpha,beta))[0]

				for j in range(0,i):

					if abs(q-1.) <= 0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

						intens[i] += alpha*np.exp(-1.*(seq[i] - seq[j]))

					#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
					elif (q > 1.):

						intens[i] += max(alpha*np.power((1+(q-1)*(seq[i]-seq[j]))+epsilon,1/(1-q)),0.)

					elif (q < 1.) and (seq[i]-seq[j]>1/(1-q)):

						intens[i] += epsilon

					elif (q < 1.) and (seq[i]-seq[j] < 1/(1-q)):

						intens[i] += max(alpha*np.power((1+(q-1)*(seq[i]-seq[j]))+epsilon,1/(1-q)),0.)				

					else:

						intens[i] += epsilon

				#print('intens: ' + repr(intens))

				intens[intens < 0.] = 0.

		print ('Loglikelihood Train GD: ' + repr(np.sum(np.nan_to_num(np.log(intens))) - compens) + '\n')

		return - np.sum(np.nan_to_num(np.log(intens))) + max(compens,0.)

	par = minimize(logGD_QEXP, [mu_0, alpha_0, q_0], method='Nelder-Mead', tol=1e-2, options={'maxiter':10})

	print('Final Parameters: '+ repr(par.x)+'\n')

	fin_llh = logGD_QEXP(par.x)

	fin_llh = (-1)*fin_llh

	QEXP_coeffs = par.x

	alpha = QEXP_coeffs[1];

	q = QEXP_coeffs[2]

	mu = QEXP_coeffs[0]

	llh_renorm_alpha = fin_llh

	llh_renorm_sqrt = fin_llh

	llh_renorm_q = fin_llh

	if abs(q-1.) <= 0.001 : #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

		statcriter = abs(alpha)

		if (statcriter >= 1.) or np.isinf(fin_llh) or np.isnan(fin_llh) or (fin_llh > 0.):

			par_renorm_alpha = [Delta*(1.-1./(1.+eps)),par.x[1]/(statcriter*(1+eps)),par.x[2]]

			par_renorm_sqrt = par_renorm_alpha

			par_renorm_q = par_renorm_alpha

			llh_renorm_alpha = logGD_QEXP(par_renorm_alpha)

			llh_renorm_sqrt = logGD_QEXP(par_renorm_sqrt)

			llh_renorm_q = logGD_QEXP(par_renorm_q)

			llh_renorm_alpha *= -1

			llh_renorm_sqrt *= -1

			llh_renorm_q *= -1			

	#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
	elif (q != 1.) and (q < 2.):

		statcriter = abs(alpha/(2-q))

		if (statcriter >= 1.) or np.isinf(fin_llh) or np.isnan(fin_llh) or (fin_llh > 0.):

			par_renorm_q = [Delta*(1.-1./(1.+eps)),par.x[1],2-(2-par.x[2])*(statcriter*(1+eps))]

			par_renorm_alpha = [Delta*(1.-1./(1.+eps)),par.x[1]/(statcriter*(1+eps)),par.x[2]]

			par_renorm_sqrt = [Delta*(1.-1./(1.+eps)),par.x[1]/np.sqrt(statcriter*(1+eps)),2-(2-par.x[2])*np.sqrt(statcriter*(1+eps))]

			llh_renorm_q = logGD_QEXP(par_renorm_q)

			llh_renorm_alpha = logGD_QEXP(par_renorm_alpha)

			llh_renorm_sqrt = logGD_QEXP(par_renorm_sqrt)

			llh_renorm_alpha *= -1

			llh_renorm_sqrt *= -1

			llh_renorm_q *= -1	

			print('par_renorm_alpha: '+repr(par_renorm_alpha))

			print('par_renorm_sqrt: ' + repr(par_renorm_sqrt))

			print('par_renorm_q: '+repr(par_renorm_q))

	else:

		q = 2/(1+eps)

		statcriter = abs(alpha/(2-q))

		if (statcriter >= 1.) or np.isinf(fin_llh) or np.isnan(fin_llh) or (fin_llh > 0.):

			par_renorm_q = [Delta*(1.-1./(1.+eps)),par.x[1],2-(2-par.x[2])*(statcriter*(1+eps))]

			par_renorm_alpha = [Delta*(1.-1./(1.+eps)),par.x[1]/(statcriter*(1+eps)),par.x[2]]

			par_renorm_sqrt = [Delta*(1.-1./(1.+eps)),par.x[1]/np.sqrt(statcriter*(1+eps)),2-(2-par.x[2])*np.sqrt(statcriter*(1+eps))]

			llh_renorm_q = logGD_QEXP(par_renorm_q)

			llh_renorm_alpha = logGD_QEXP(par_renorm_alpha)

			llh_renorm_sqrt = logGD_QEXP(par_renorm_sqrt)

			llh_renorm_alpha *= -1

			llh_renorm_sqrt *= -1

			llh_renorm_q *= -1	

			print('par_renorm_alpha: '+repr(par_renorm_alpha))

                        print('par_renorm_sqrt: ' + repr(par_renorm_sqrt))

                        print('par_renorm_q: '+repr(par_renorm_q))


	print('QEXP_statcriter: ' + repr(statcriter))

	print('llh_renorm_alpha:' + repr(llh_renorm_alpha))

	print('llh_renorm_sqrt: '+repr(llh_renorm_sqrt))

	print('llh_renorm_q: ' + repr(llh_renorm_q))

	K1_Param = {'QEXP_coeffs': par.x, 'K1_Type': 'QEXP', 'QEXP_statcriter': statcriter, 'final_llh': fin_llh, 'llh_renorm_alpha': llh_renorm_alpha, \
	'llh_renorm_q': llh_renorm_q, 'llh_renorm_sqrt': llh_renorm_sqrt}

	return K1_Param
