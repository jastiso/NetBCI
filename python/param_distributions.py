#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 13:17:03 2018

Plot parameter histograms

@author: stiso
"""

import numpy as np
from pathlib import Path
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import plotly.plotly as py

top_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/gc/'
subjs = ['001',  '002','003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018','019', '020']
bands = ['alpha', 'beta', 'low_gamma', 'gamma', 'low_gamma']
sens = ['grad', 'mag']



for e in sens:
    alphas = []
    betas = []
    ranks = []
    alphas_bl = []
    betas_bl = []
    ranks_bl = []
    for s in subjs:
        for b in bands: 
            # Load
            my_file = Path("".join([top_dir, s, '/', e, '/gc_', b ,'_params.npy']))
            if my_file.is_file():
                params = np.load("".join([top_dir, s, '/', e, '/gc_', b ,'_params.npy'])).item()
                print(params['alpha']) 
                alphas.append(params["alpha"])
                print(params['beta'])
                betas.append(params["beta"])
                print(params['rank'])
                ranks.append(params["rank"])
                               
        
#Now for Baseline
#    for s in subjs:
#        for b in bands:
#            my_file = Path("".join([top_dir, s, '/', e, '/ac_', b ,'_params_bl.npy']))
#            if my_file.is_file():
#                # Load
#                params = np.load("".join([top_dir, s, '/', e, '/ac_', b ,'_params_bl.npy'])).item()
#                print(params['alpha']) 
#                alphas_bl.append(params["alpha"])
#                print(params['beta'])
#                betas_bl.append(params["beta"])
#                print(params['rank'])
#                ranks_bl.append(params["rank"])
                
                
                # Plot
                
    fig = plt.figure()

    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(alphas, num_bins, normed=1, facecolor='green', alpha=0.5)
    plt.xlabel('alpha')
    plt.ylabel('Frequency')
    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/gc/', e, 'alpha_hist.png']))

    fig = plt.figure()

    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(betas, num_bins, normed=1, facecolor='blue', alpha=0.5)
    plt.xlabel('beta')
    plt.ylabel('Frequency')
    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/gc/', e, 'beta_hist.png']))

    fig = plt.figure()

    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(ranks, num_bins, normed=1, facecolor='purple', alpha=0.5)
    plt.xlabel('rank')
    plt.ylabel('Frequency')
    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/gc/', e, 'rank_hist.png']))

    # Plot baseline
                    
#    fig = plt.figure()

#    num_bins = 50
#    # the histogram of the data
#    n, bins, patches = plt.hist(alphas_bl, num_bins, normed=1, facecolor='grey', alpha=0.5)
#    plt.xlabel('alpha')
#    plt.ylabel('Frequency')
#    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/ac/', e, 'alpha_hist_bl.png']))

#    fig = plt.figure()

#    num_bins = 50
    # the histogram of the data
#    n, bins, patches = plt.hist(betas_bl, num_bins, normed=1, facecolor='grey', alpha=0.5)
#    plt.xlabel('beta')
#    plt.ylabel('Frequency')
#    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/ac/', e, 'beta_hist_bl.png']))

 #   fig = plt.figure()

#    num_bins = 50
#    # the histogram of the data
 #   n, bins, patches = plt.hist(ranks_bl, num_bins, normed=1, facecolor='grey', alpha=0.5)
 #   plt.xlabel('rank')
 #   plt.ylabel('Frequency')
 #   fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/ac/', e, 'rank_hist_bl.png']))

