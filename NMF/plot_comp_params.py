#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 13:21:34 2018

@author: stiso
"""

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

top_dir1 = '/Users/stiso/Documents/Python/NetBCI/NMF/param/'
top_dir2 = '/Users/stiso/Documents/Python/NetBCI/NMF/'
subjs = ['001',  '002','003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018','019', '020']
bands = ['alpha', 'beta', 'low_gamma']
sens = ['grad']



for e in sens:
    for b in bands:
        print(b)
        alphas1 = []
        betas1 = []
        ranks1 = []
        alphas2 = []
        betas2 = []
        ranks2 = []
        errors1 = []
        errors2 = []
        for s in subjs: 
            # Load
            params1 = np.load("".join([top_dir1, s, '/', e, 'wpli_', b ,'_params.npy'])).item()
            alphas1.append(params1["alpha"])
            betas1.append(params1["beta"])
            ranks1.append(params1["rank"])
            params2 = np.load("".join([top_dir2, s, '/', e, '/wpli_', b ,'_params.npy'])).item()
            alphas2.append(params2["alpha"])
            betas2.append(params2["beta"])
            ranks2.append(params2["rank"])
            err1 = np.load("".join([top_dir1, s, '/', e, 'wpli_', b ,'_err.npy']))
            err2 = np.load("".join([top_dir2, s, '/', e, '/wpli_', b ,'_err.npy']))
            errors1.append(err1[err1.size-1])
            errors2.append(err2[err2.size-1])

        # plot
        fig = plt.figure()
        plt.scatter(alphas1, alphas2, facecolor = 'green')
        plt.xlabel('param1')
        plt.ylabel('param2')
        # y=x line
        plt.plot([np.min(alphas1 + alphas2),np.max(max(alphas1 + alphas2))], [np.min(alphas1 + alphas2),np.max(max(alphas1 + alphas2))], color = "black")
        fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/', b, '_alpha_wpli.png']))
        
        fig = plt.figure()
        plt.scatter(betas1, betas2, facecolor = 'blue')
        plt.xlabel('param1')
        plt.ylabel('param2')
        # y=x line
        plt.plot([np.min(betas1 + betas2),np.max(max(betas1 + betas2))], [np.min(betas1 + betas2),np.max(max(betas1 + betas2))], color = "black")
        fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/', b, '_beta_wpli.png']))
         
        fig = plt.figure()
        plt.scatter(ranks1, ranks2, facecolor = 'purple')
        plt.xlabel('param1')
        plt.ylabel('param2')
        # y=x line
        plt.plot([np.min(ranks1 + ranks2),np.max(max(ranks1 + ranks2))], [np.min(ranks1 + ranks2),np.max(max(ranks1 + ranks2))], color = "black")
        fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/', b, '_rank_wpli.png']))
        
        fig = plt.figure()
        plt.scatter(errors1, errors2, facecolor = 'orange')
        plt.xlabel('param1')
        plt.ylabel('param2')
        # y=x line
        plt.plot([np.min(errors1 + errors2),np.max(max(errors1 + errors2))], [np.min(errors1 + errors2),np.max(max(errors1 + errors2))], color = "black")
        fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/', b, '_err_wpli.png']))

