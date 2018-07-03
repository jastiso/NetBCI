#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 11:34:50 2018

Run NMF on a single subject
(i) run the NMF algorithm R times per network configu-ration matrix, (ii) concatenate subgraph matrixWacross R runs
into an aggregate matrix with dimensions E × (R ∗ ¯
m), and (iii)
into an aggregate matrix with dimensions E × (R ∗ ¯
m), and (iii)
apply NMF to the aggregate matrix to determine a final set of
subgraphs and expression coefficients.

@author: stiso
"""
import os
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
os.chdir('/Users/stiso/Documents/Python/Echobase-master/')
#from Echobase import optimize_nmf
import optimize_nmf

#%% Global variables

os.chdir('/Users/stiso/Documents/Python/NetBCI/')
subj = '003';
band = 'alpha'
eType = 'mag'

data_dir = '/Users/stiso/Documents/Matlab/NetBCI/NMF/'
save_dir = ''.join(['/Users/stiso/Documents/Python/NetBCI/NMF/', subj, '/'])
# make directory
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
#%% Load data 

data = io.loadmat(''.join([data_dir, 'gc_', band, '_', eType, '_', subj, '.mat']))
# should be nWin x nCon
val = np.transpose(np.array(data['A']))
nWin = val.shape[0]

#%% Test all parameters

alpha_range = (.001,2)
beta_range = (.001,2)
m_range = (2,20)
n_param = 10000

#%% get fold id
k = 10
fold_id = [[] for i in range(k)]
indices = np.arange(0,nWin-1)
# get sizes of each fold
sizes = np.tile(int(np.floor(nWin/k)), k)
for i in range(nWin%k):
    sizes[i] = sizes[i] + 1
# make ids
indices = np.random.choice(nWin, replace = False, size=nWin)
for i in range(k):
    fold_id[i] = indices[0:sizes[i]].tolist()
    indices = np.delete(indices, np.arange(0,sizes[i]))
    
params = optimize_nmf.gen_random_sampling_paramset(m_range, alpha_range, beta_range, n_param, fold_id)

#%% Cross validation for parameters

qual_meas = [];
for i in range(n_param):
    qual_meas.append(optimize_nmf.run_xval_paramset(val, params[i]))

#%% Use the best parameters

[opt_dict, opt_param] = optimize_nmf.find_optimum_xval_paramset(params, qual_meas)

#%% Consensus clustering for best parameters

n_proc = 1
opt_alpha = opt_param.get('alpha')
opt_beta = opt_param.get('beta')
opt_rank= opt_param.get('rank')
n_seed = 100;
[subset, coeff, err] = optimize_nmf.consensus_nmf(val, opt_alpha, opt_beta,
                  opt_rank, n_seed, n_proc,)
np.save("".join([save_dir, band, '_subset']), subset)
np.save("".join([save_dir, band, '_coeff']), coeff)
np.save("".join([save_dir, band, '_err']), err)

#%% Plot

#plt.imshow(coeff, cmap='plasma',  ax=plt.subplots(figsize=(10,10)))
#plt.show()
fig = plt.figure()
plt.plot(np.transpose(coeff))
fig.savefig("".join([save_dir, band, '_', subj, '_node_exp.png']))

fig = plt.figure()
plt.plot(coeff)
fig.savefig("".join([save_dir, band, '_', subj, '_sg_exp.png']))

coeff_nb = coeff[0:6,0:20706]
fig = plt.figure()
plt.plot(np.transpose(coeff_nb))
fig.savefig("".join([save_dir, band, '_', subj, '_node_exp_nb.png']))

fig = plt.figure()
plt.plot(coeff_nb)
fig.savefig("".join([save_dir, band, '_', subj, '_sg_exp_nb.png']))
