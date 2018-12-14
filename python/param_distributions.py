#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 13:17:03 2018

Plot parameter histograms

@author: stiso
"""

import numpy as np
top_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/'
subjs = ['003', '004', '005']
bands = ['theta', 'alpha', 'beta', 'low_gamma', 'gamma']
sens = ['mag']

for s in subjs:
    for b in bands: 
        # Load
        params = np.load("".join([top_dir, s, '/', b ,'_params.npy'])).item()
        print(params['alpha']) 
        print(params['beta'])
        
#Now for Baseline
for s in subjs:
    for b in bands: 
        # Load
        params = np.load("".join([top_dir, s, '/', b ,'_params_bl.npy'])).item()
        print(params['alpha']) 
        print(params['beta'])