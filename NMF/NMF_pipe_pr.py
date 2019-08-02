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
import sys
#os.chdir('/Users/stiso/Documents/Python/Echobase-master/')
#from Echobase import optimize_nmf
import optimize_nmf
import importlib.util
#spec = importlib.util.spec_from_file_location('optimize_nmf', '/data/jag/bassett-lab/jstiso/Echobase-master/Echobase/Network/Partitioning/Subgraph/optimize_nmf.py')
#optimize_nmf = importlib.util.module_from_spec(spec)
#spec.loader.exec_module(optimize_nmf)

def pipe(subj, band, eType):

#%% Global variables
    os.chdir('/data/jag/bassett-lab/jstiso/Python/NetBCI/')
    data_dir = '/data/jag/bassett-lab/jstiso/Python/NetBCI/data/'
    save_dir = ''.join(['/data/jag/bassett-lab/jstiso/Python/NetBCI/NMF/', subj, '/', eType, '/'])
# make directory
    if not os.path.exists(save_dir):
        os.makedirs(save_dir, exist_ok=True)
#%% Load data 

    data = io.loadmat(''.join([data_dir, 'pr_gc_', band, '_', eType, '_', subj, '.mat']))
    val = np.transpose(np.array(data['A']))
    nWin = val.shape[0]

#%% Test all parameters

    alpha_range = (0.01,1.5)
    beta_range = (0.01,1.5)
    m_range = (2,15)
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
    np.save("".join([save_dir,'pr_gc_', band, '_params']), opt_param)
    
#%% Consensus clustering for best parameters

    n_proc = 1
    opt_alpha = opt_param.get('alpha')
    opt_beta = opt_param.get('beta')
    opt_rank= opt_param.get('rank')
    n_seed = 100;
    [subset, coeff, err] = optimize_nmf.consensus_nmf(val, opt_alpha, opt_beta, opt_rank, n_seed, n_proc)
    np.save("".join([save_dir, 'pr_gc_',band, '_subset']), subset)
    np.save("".join([save_dir,'pr_gc_', band, '_coeff']), coeff)
    np.save("".join([save_dir, 'pr_gc_',band, '_err']), err)

#%% work with command line

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    pipe(*sys.argv[1:])



