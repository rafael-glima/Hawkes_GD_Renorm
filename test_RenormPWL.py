
import scipy.io
from scipy.optimize import minimize
import numpy as np
from scipy.integrate import quad
from scipy.special import lambertw

class par():

	x = [5.,1.1,1.1,1.1]

eps = 0.1

PWL_statcriter = par.x[0]/((par.x[2]-1)*np.power(par.x[1],par.x[2]-1))

par_renorm_K = [par.x[0]/(PWL_statcriter*(1+eps)),par.x[1],par.x[2],par.x[3]]

par_renorm_c = [par.x[0],np.power((1+eps)*PWL_statcriter,1/par.x[2])*par.x[1],par.x[2],par.x[3]]

Delta_p = PWL_statcriter*(1+eps)*(par.x[2]-1)*np.power(par.x[1],par.x[2]-1.)

print('Delta_p: '+repr(Delta_p))

par_renorm_p = [par.x[0],par.x[1],1+ lambertw(Delta_p*np.log(par.x[1])).real/np.log(par.x[1]),par.x[3]]

par_renorm_Kc = [par.x[0]/np.sqrt(PWL_statcriter*(1+eps)),par.x[1]*np.power(PWL_statcriter*(1+eps),1/(2*(par.x[2]-1))),par.x[2],par.x[3]]

Delta_Kp = np.sqrt(PWL_statcriter*(1+eps))*(par.x[2]-1)*np.power(par.x[1],par.x[2]-1.)

par_renorm_Kp = [par.x[0]/np.sqrt(PWL_statcriter*(1+eps)),par.x[1],1+lambertw(Delta_Kp*np.log(par.x[1])).real/np.log(par.x[1]),par.x[3]]

print('par_renorm_K:' + repr(par_renorm_K))

print('par_renorm_c:'+repr(par_renorm_c))

print('par_renorm_p:'+repr(par_renorm_p))

print('par_renorm_Kc:'+repr(par_renorm_Kc))

print('par_renorm_Kp:'+repr(par_renorm_Kp))
