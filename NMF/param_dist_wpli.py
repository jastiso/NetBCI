#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 10:36:19 2018

@author: stiso
"""


import numpy as np
from pathlib import Path
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import plotly.plotly as py

top_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/'
subjs = ['001',  '002','003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018','019', '020']
bands = ['alpha', 'beta', 'low_gamma','low_gamma']
sens = ['grad']



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
            my_file = Path("".join([top_dir, s, '/', e, 'wpli_', b ,'_params.npy']))
            if my_file.is_file():
                params = np.load("".join([top_dir, s, '/', e, 'wpli_', b ,'_params.npy'])).item()
                print(params['alpha']) 
                alphas.append(params["alpha"])
                print(params['beta'])
                betas.append(params["beta"])
                print(params['rank'])
                ranks.append(params["rank"])
                               
                
                # Plot
                
    fig = plt.figure()

    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(alphas, num_bins, normed=1, facecolor='green', alpha=0.5)
    plt.xlabel('alpha')
    plt.ylabel('Frequency')
    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/', e, 'alpha_hist_wpli.png']))

    fig = plt.figure()

    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(betas, num_bins, normed=1, facecolor='blue', alpha=0.5)
    plt.xlabel('beta')
    plt.ylabel('Frequency')
    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/', e, 'beta_hist_wpli.png']))

    fig = plt.figure()

    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(ranks, num_bins, normed=1, facecolor='purple', alpha=0.5)
    plt.xlabel('rank')
    plt.ylabel('Frequency')
    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/', e, 'rank_hist_wpli.png']))

print(np.mean(ranks))
print(np.mean(betas))
print(np.mean(alphas))