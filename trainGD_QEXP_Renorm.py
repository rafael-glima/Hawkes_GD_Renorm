import scipy.io
from scipy.optimize import minimize
import numpy as np
from scipy.integrate import quad
import math
#import numpy.random as np.random

def trainGD_QEXP(seq,eps):

	alpha_0 = np.random.rand()
	beta_0 = np.random.rand() #2*alpha_0
	q_0 = 1. + np.random.rand()
	mu_0 = np.random.rand()

	# input_data = scipy.io.loadmat('4Kern_newSNS_10seqT100000_highMuandfreq0.15.mat')
	# seq = input_data['Seq2'][1][0][0]

	# seq = seq[:300]

	T = seq[-1]#-seq[0]
	Delta = len(seq)/T

	#seq = np.array([1.,2.,3.,4.,5.])

	bnds = ((0,None),(0,None),(0,None))

	def logGD_QEXP(QEXP_coeffs):

		def funcqexp(x,alpha,beta,q):

			print('q: '+ repr(q))

			if abs(q-1.)<=0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

				return alpha*np.exp(-beta*x)

			elif (q != 1.) and (1 + (1-q)*beta*x > 0):

				return np.power(alpha*(1+(q-1)*beta*x),1/(1-q))

			else:

				return 0

		alpha = QEXP_coeffs[1];

		beta = QEXP_coeffs[2];

		q = QEXP_coeffs[3]

		mu = QEXP_coeffs[0]


		if abs(q-1.) <= 0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

			if (alpha/beta < 1.) and (alpha/beta >= 0.):

				mu = mu #(1.-alpha/beta)*Delta;

			else:

				#mu = 0. 
				return np.inf

		#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
		elif (q != 1.) and ((q-1)*beta > 0.):

			if (alpha*(q-1)/(beta*(2-q)) < 1.) and (alpha*(q-1)/(beta*(2-q)) > 0.):

				mu = mu#(1.- alpha*(q-1)/(beta*(2-q)))*Delta;

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

				compens += (alpha/beta)*(1-np.exp(-beta*(T-seq[i])))

			#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
			elif (q != 1.) and ((q-1)*beta > 0.):

				compens += alpha*((1-q)/(beta*(2-q)))*(np.power(1+(q-1)*beta*(T-seq[i]), (2-q)/(1-q))-1)

			else:

				compens += 0.

			#compens += (alpha/beta)*(1-np.exp(-beta*(T-seq[i])))#quad(funcexp,0,T-seq[i], args=(alpha,beta))[0]

			for j in range(0,i):

				if abs(q-1.) <= 0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

					intens[i] += alpha*np.exp(-beta*(seq[i] - seq[j]))

				#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
				elif (q != 1.) and ((q-1)*beta > 0.):

					intens[i] += alpha*((1-q)/(beta*(2-q)))*np.power((1+(q-1)*beta*(seq[i]-seq[j])),(2-q)/(1-q))

				else:

					intens[i] += 0.

		intens[intens < 0.] = 0.	

		print ('Loglikelihood Train GD: ' + repr(np.sum(np.nan_to_num(np.log(intens))) - compens) + '\n')

		return - np.sum(np.nan_to_num(np.log(intens))) + compens

	par = minimize(logGD_QEXP, [mu_0, alpha_0, beta_0, q_0], method='Nelder-Mead', tol=1e-2, options={'maxiter':10})

	print('Final Parameters: '+ repr(par.x)+'\n')

	fin_llh = logGD_QEXP(par.x)

	fin_llh = (-1)*fin_llh

	QEXP_coeffs = par.x

	alpha = QEXP_coeffs[1];

	beta = QEXP_coeffs[2];

	q = QEXP_coeffs[3]

	mu = QEXP_coeffs[0]

	llh_renorm_alpha = fin_llh

	llh_renorm_beta = fin_llh

	llh_renorm_sqrt = fin_llh

	llh_renorm_q = fin_llh

	if abs(q-1.) <= 0.001: #math.isclose(q, 1., rel_tol=1e-3, abs_tol=0.0):

		statcriter = alpha/beta

		if statcriter >= 1.:

			par_renorm_alpha = [Delta*(1-1/(1.+eps)),par.x[1]/(1+eps),par.x[2],par.x[3]]

			par_renorm_beta = [Delta*(1-1/(1.+eps)),par.x[1],par.x[2]*(1+eps),par.x[3]]

			par_renorm_sqrt = [Delta*(1-1/(1.+eps)),par.x[1]/np.sqrt(1+eps),par.x[2]*np.sqrt(1+eps),par.x[3]]

			par_renorm_q = par_renorm_sqrt

			llh_renorm_alpha = llhGD_QEXP(par_renorm_alpha)

			llh_renorm_beta = llhGD_QEXP(par_renorm_beta)

			llh_renorm_sqrt = llhGD_QEXP(par_renorm_sqrt)

			llh_renorm_q = llh_renorm_sqrt

			llh_renorm_alpha *= -1

			llh_renorm_beta *= -1

			llh_renorm_sqrt *= -1

			llh_renorm_q *= -1			

	#elif (q != 1.) and (1 + (q-1)*beta*x > 0.):
	elif (q != 1.) and ((q-1)*beta > 0.):

		statcriter = (q-1)/(2-q)

		if statcriter >= 1.:

			par_renorm_q = [Delta*(1-1/(1+eps)),par.x[1],par.x[2],(par.x[3]*(statcriter+1)+2+eps-2*statcriter)/(2+eps)]

			par_renorm_alpha = par_renorm_q

			par_renorm_beta = par_renorm_q

			par_renorm_sqrt = par_renorm_q

			llh_renorm_q = llhGD_QEXP(par_renorm_q)

			llh_renorm_alpha = llh_renorm_q

			llh_renorm_beta = llh_renorm_q

			llh_renorm_sqrt = llh_renorm_q

			llh_renorm_alpha *= -1

			llh_renorm_beta *= -1

			llh_renorm_sqrt *= -1

			llh_renorm_q *= -1	

	else:

		statcriter = 0.
	print('QEXP_statcriter: ' + repr(statcriter))

	K1_Param = {'QEXP_coeffs': par.x, 'K1_Type': 'QEXP', 'QEXP_statcriter': statcriter, 'final_llh': fin_llh, 'llh_renorm_alpha': llh_renorm_alpha, \
	'llh_renorm_beta': llh_renorm_beta,'llh_renorm_q': llh_renorm_q, 'llh_renorm_sqrt': llh_renorm_sqrt}

	return K1_Param
