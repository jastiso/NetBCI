
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


def pipe():
#%% Global variables

    data_dir = '/data/jag/bassett-lab/jstiso/Python/NetBCI/data/'
    save_dir = '/data/jag/bassett-lab/jstiso/Python/NetBCI/NMF/'

#%% Load data 

    rand_data = io.loadmat('/data/jag/bassett-lab/jstiso/Python/NetBCI/data/rand_A.mat')
    val_rand = np.transpose(np.array(rand_data['A']))
    nWin = val_rand.shape[0]
    reg_data = io.loadmat(''.join([data_dir, 'reg_A.mat']))
    val_reg = np.transpose(np.array(reg_data['A']))


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

    qual_meas_rand = [];
    qual_meas_reg = [];
    for i in range(n_param):
        qual_meas_rand.append(optimize_nmf.run_xval_paramset(val_rand, params[i]))
        #qual_meas_reg.append(optimize_nmf.run_xval_paramset(val_reg, params[i]))

#%% Use the best parameters

    [opt_dict_rand, opt_param_rand] = optimize_nmf.find_optimum_xval_paramset(params, qual_meas_rand)
    np.save("".join([save_dir,'rand_params']), opt_param_rand)
    #[opt_dict_reg, opt_param_reg] = optimize_nmf.find_optimum_xval_paramset(params, qual_meas_reg)
    #np.save("".join([save_dir,'reg_params']), opt_param_reg)  

#%% Consensus clustering for best parameters

    n_proc = 1
    opt_alpha = opt_param_rand.get('alpha')
    opt_beta = opt_param_rand.get('beta')
    opt_rank= opt_param_rand.get('rank')
    n_seed = 100;
    [subset, coeff, err] = optimize_nmf.consensus_nmf(val_rand, opt_alpha, opt_beta, opt_rank, n_seed, n_proc)
    np.save("".join([save_dir, 'rand_subset']), subset)
    np.save("".join([save_dir,'rand_coeff']), coeff)
    np.save("".join([save_dir, 'rand_err']), err)

    #opt_alpha = opt_param_reg.get('alpha')
    #opt_beta = opt_param_reg.get('beta')
    #opt_rank= opt_param_reg.get('rank')
    #[subset, coeff, err] = optimize_nmf.consensus_nmf(val_reg, opt_alpha, opt_beta, opt_rank, n_seed, n_proc)
    #np.save("".join([save_dir, 'reg_subset']), subset)
    #np.save("".join([save_dir,'reg_coeff']), coeff)
    #np.save("".join([save_dir, 'reg_err']), err)

#%% work with command line

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    pipe()



