
import numpy as np
import scipy.io
import math
import matplotlib.pyplot as plt
import random as rd
rd.seed(2018)
np.random.seed(2018)

import pandas as pd
import datetime as DT
import time

from joblib import Parallel, delayed
import multiprocessing
import parmap
from multiprocessing import Pool, freeze_support
import itertools

from trainGD_EXP_Renorm import trainGD_EXP
from trainGD_PWL_Renorm import trainGD_PWL
from trainGD_QEXP_Renorm_New import trainGD_QEXP
from trainGD_RAY_Renorm import trainGD_RAY

input_data = scipy.io.loadmat('4Kern_Renorm_10seq_T30000.mat')

eps = 0.001

n_of_seq = 10

resolution = 40

llh_GD_EXP = np.array([])
llh_GD_PWL = np.array([])
llh_GD_QEXP = np.array([])
llh_GD_RAY = np.array([])

llh_GD_EXP_Renorm_alpha = np.array([])
llh_GD_PWL_Renorm_K = np.array([])
llh_GD_QEXP_Renorm_a = np.array([])
llh_GD_RAY_Renorm_gamma = np.array([])

llh_GD_EXP_Renorm_beta = np.array([])
llh_GD_PWL_Renorm_c = np.array([])
llh_GD_PWL_Renorm_p = np.array([])
llh_GD_QEXP_Renorm_q = np.array([])
llh_GD_RAY_Renorm_eta = np.array([])

llh_GD_EXP_Renorm_alphabeta = np.array([])
llh_GD_PWL_Renorm_Kc = np.array([])
llh_GD_QEXP_Renorm_aq = np.array([])
llh_GD_RAY_Renorm_gammaeta = np.array([])

llh_GD_PWL_Renorm_Kp = np.array([])

for i in range(1,5):

	ind_seq = 'Seq' + repr(i)

	print(ind_seq)

	for j in range(0,n_of_seq):

		print(j)

		seq = input_data[ind_seq][j][0][0]

		EXP_Param = trainGD_EXP(seq,eps)

		PWL_Param = trainGD_PWL(seq,eps)

		QEXP_Param = trainGD_QEXP(seq,eps)

		RAY_Param = trainGD_RAY(seq,eps)

		llh_GD_EXP = np.append(llh_GD_EXP,EXP_Param['final_llh'])
		llh_GD_PWL = np.append(llh_GD_PWL,PWL_Param['final_llh'])
		llh_GD_QEXP = np.append(llh_GD_QEXP,QEXP_Param['final_llh'])
		llh_GD_RAY = np.append(llh_GD_RAY,RAY_Param['final_llh'])

		llh_GD_EXP_Renorm_alpha = np.append(llh_GD_EXP_Renorm_alpha,EXP_Param['llh_renorm_alpha'])
		llh_GD_PWL_Renorm_K = np.append(llh_GD_PWL_Renorm_K,PWL_Param['llh_renorm_K'])
		llh_GD_QEXP_Renorm_a = np.append(llh_GD_QEXP_Renorm_a,QEXP_Param['llh_renorm_alpha'])
		llh_GD_RAY_Renorm_gamma = np.append(llh_GD_RAY_Renorm_gamma,RAY_Param['llh_renorm_gamma'])

		llh_GD_EXP_Renorm_beta = np.append(llh_GD_EXP_Renorm_beta,EXP_Param['llh_renorm_beta'])
		llh_GD_PWL_Renorm_c = np.append(llh_GD_PWL_Renorm_c,PWL_Param['llh_renorm_c'])
		llh_GD_PWL_Renorm_p = np.append(llh_GD_PWL_Renorm_p,PWL_Param['llh_renorm_p'])
		llh_GD_QEXP_Renorm_q = np.append(llh_GD_QEXP_Renorm_q,QEXP_Param['llh_renorm_q'])
		llh_GD_RAY_Renorm_eta = np.append(llh_GD_RAY_Renorm_eta,RAY_Param['llh_renorm_eta'])

		llh_GD_EXP_Renorm_alphabeta = np.append(llh_GD_EXP_Renorm_alphabeta,EXP_Param['llh_renorm_sqrt'])
		llh_GD_PWL_Renorm_Kc = np.append(llh_GD_PWL_Renorm_p,PWL_Param['llh_renorm_Kc'])
		llh_GD_QEXP_Renorm_aq = np.append(llh_GD_QEXP_Renorm_aq,QEXP_Param['llh_renorm_sqrt'])
		llh_GD_RAY_Renorm_gammaeta = np.append(llh_GD_RAY_Renorm_gammaeta,RAY_Param['llh_renorm_sqrt'])

		llh_GD_PWL_Renorm_Kp = np.append(llh_GD_PWL_Renorm_Kp,PWL_Param['llh_renorm_Kp'])


		
print('llh_GD_EXP: ' + repr(llh_GD_EXP) + '\n')
print('llh_GD_EXP_Renorm_alpha: ' + repr(llh_GD_EXP_Renorm_alpha) + '\n')
print('llh_GD_EXP_Renorm_beta: ' + repr(llh_GD_EXP_Renorm_beta) + '\n')
print('llh_GD_EXP_Renorm_alphabeta: ' + repr(llh_GD_EXP_Renorm_alphabeta) + '\n')

print('llh_GD_PWL: ' + repr(llh_GD_PWL) + '\n')
print('llh_GD_PWL_Renorm_K: ' + repr(llh_GD_PWL_Renorm_K) + '\n')
print('llh_GD_PWL_Renorm_c: ' + repr(llh_GD_PWL_Renorm_c) + '\n')
print('lllh_GD_PWL_Renorm_p: ' + repr(llh_GD_PWL_Renorm_p) + '\n')
print('llh_GD_PWL_Renorm_Kc: ' + repr(llh_GD_PWL_Renorm_Kc) + '\n')
print('llh_GD_PWL_Renorm_Kp: ' + repr(llh_GD_PWL_Renorm_Kp) + '\n')

print('llh_GD_QEXP: ' + repr(llh_GD_QEXP) + '\n')
print('llh_GD_QEXP_Renorm_a: ' + repr(llh_GD_QEXP_Renorm_a) + '\n')
print('lllh_GD_QEXP_Renorm_q: ' + repr(llh_GD_QEXP_Renorm_q) + '\n')
print('llh_GD_QEXP_Renorm_aq: ' + repr(llh_GD_QEXP_Renorm_aq) + '\n')

print('llh_GD_RAY: ' + repr(llh_GD_RAY) + '\n')
print('llh_GD_RAY_Renorm_gamma: ' + repr(llh_GD_RAY_Renorm_gamma) + '\n')
print('llh_GD_RAY_Renorm_eta: ' + repr(llh_GD_RAY_Renorm_eta) + '\n')
print('llh_GD_RAY_Renorm_gammaeta: ' + repr(llh_GD_RAY_Renorm_gammaeta) + '\n')

f = open('Exp_Synthetic_Renorm_T30000_eps_0.001.txt','w')
f.write('llh_GD_EXP: ' + repr(llh_GD_EXP) + '\n')
f.write('llh_GD_EXP_Renorm_alpha: ' + repr(llh_GD_EXP_Renorm_alpha) + '\n')
f.write('llh_GD_EXP_Renorm_beta: ' + repr(llh_GD_EXP_Renorm_beta) + '\n')
f.write('llh_GD_EXP_Renorm_alphabeta: ' + repr(llh_GD_EXP_Renorm_alphabeta) + '\n')

f.write('llh_GD_PWL: ' + repr(llh_GD_PWL) + '\n')
f.write('llh_GD_PWL_Renorm_K: ' + repr(llh_GD_PWL_Renorm_K) + '\n')
f.write('llh_GD_PWL_Renorm_c: ' + repr(llh_GD_PWL_Renorm_c) + '\n')
f.write('lllh_GD_PWL_Renorm_p: ' + repr(llh_GD_PWL_Renorm_p) + '\n')
f.write('llh_GD_PWL_Renorm_Kc: ' + repr(llh_GD_PWL_Renorm_Kc) + '\n')
f.write('llh_GD_PWL_Renorm_Kp: ' + repr(llh_GD_PWL_Renorm_Kp) + '\n')

f.write('llh_GD_QEXP: ' + repr(llh_GD_QEXP) + '\n')
f.write('llh_GD_QEXP_Renorm_a: ' + repr(llh_GD_QEXP_Renorm_a) + '\n')
f.write('lllh_GD_QEXP_Renorm_q: ' + repr(llh_GD_QEXP_Renorm_q) + '\n')
f.write('llh_GD_QEXP_Renorm_aq: ' + repr(llh_GD_QEXP_Renorm_aq) + '\n')

f.write('llh_GD_RAY: ' + repr(llh_GD_RAY) + '\n')
f.write('llh_GD_RAY_Renorm_gamma: ' + repr(llh_GD_RAY_Renorm_gamma) + '\n')
f.write('llh_GD_RAY_Renorm_eta: ' + repr(llh_GD_RAY_Renorm_eta) + '\n')
f.write('llh_GD_RAY_Renorm_gammaeta: ' + repr(llh_GD_RAY_Renorm_gammaeta) + '\n')

f.close()

np.savez('llh_arrays_eps_0.001_T30000.npz', llh_GD_EXP=llh_GD_EXP,llh_GD_EXP_Renorm_alpha=llh_GD_EXP_Renorm_alpha,llh_GD_EXP_Renorm_beta=llh_GD_EXP_Renorm_beta,llh_GD_EXP_Renorm_alphabeta=llh_GD_EXP_Renorm_alphabeta,llh_GD_PWL=llh_GD_PWL,llh_GD_PWL_Renorm_K=llh_GD_PWL_Renorm_K,llh_GD_PWL_Renorm_c=llh_GD_PWL_Renorm_c,llh_GD_PWL_Renorm_p=llh_GD_PWL_Renorm_p,llh_GD_PWL_Renorm_Kc=llh_GD_PWL_Renorm_Kc,llh_GD_PWL_Renorm_Kp=llh_GD_PWL_Renorm_Kp,llh_GD_QEXP=llh_GD_QEXP,llh_GD_QEXP_Renorm_a=llh_GD_QEXP_Renorm_a,llh_GD_QEXP_Renorm_q=llh_GD_QEXP_Renorm_q,llh_GD_QEXP_Renorm_aq=llh_GD_QEXP_Renorm_aq,llh_GD_RAY=llh_GD_RAY,llh_GD_RAY_Renorm_gamma=llh_GD_RAY_Renorm_gamma,llh_GD_RAY_Renorm_eta=llh_GD_RAY_Renorm_eta,llh_GD_RAY_Renorm_gammaeta=llh_GD_RAY_Renorm_gammaeta)
