#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 11:08:49 2018

@author: stiso

Make plots of dynamics functional connectivity for fig 1 for net BCI project

"""

import numpy as np
from pathlib import Path
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import plotly.plotly as py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.io as sio
import seaborn as sns

top_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/'
subjs = ['002']
bands = ['beta']
sens = ['grad']
nNode = 102

#%%Plot subgraphs
for e in sens:
    for s in subjs:
        for b in bands: 
            # Load
            my_file = Path("".join([top_dir, s, '/', e, '/wpli_', b ,'_subset.npy']))
            if my_file.is_file():
                
                # now plot subgraphs
                subset = np.array(np.load("".join([top_dir, s, '/', e, '/wpli_', b ,'_subset.npy'])))
                high = np.max(subset[:,subset.shape[1]-1])
                low = np.min(subset[:,subset.shape[1]-1])
                for i in range(subset.shape[0]):
                    A = np.zeros([nNode,nNode])
                    cnt = -1
                    for col in range(nNode):
                        for row in range(col):
                            if (row != col):
                                cnt = cnt + 1;
                                A[row,col] = subset[i,cnt]
                    A = A + np.transpose(A)

                    fig = plt.figure()
                    cax = plt.imshow(A, cmap='viridis', interpolation='nearest')
                    fig.colorbar(cax)
                    fig.show()
                    if (subset[i,subset.shape[1]-1] == high):
                        fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/wpli/', e, '_', s, '_', b, '_subset_high.pdf']), format = 'pdf')
                    elif (subset[i,subset.shape[1]-1] == low):
                        fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/wpli/', e, '_', s, '_', b, '_subset_low.pdf']), format = 'pdf')
                    else:
                         fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/wpli/', e, '_', s, '_', b, '_subset_', str(i), '.pdf']), format = 'pdf') 
                         
                         
#%% Plot Connection Matrix
                         
for e in sens:
    for s in subjs:
        for b in bands: 
            # Load
            my_file = Path("".join(['/Users/stiso/Documents/MATLAB/NetBCI/NMF/wpli_', b, '_', e, '_', s, '.mat']))
            if my_file.is_file():
                
                # now plot subgraphs
                connection_mat = sio.loadmat("".join(['/Users/stiso/Documents/MATLAB/NetBCI/NMF/wpli_', b, '_', e, '_', s, '.mat']))
                A = connection_mat.get('A')

                #fig = plt.figure(figsize=(8, 6))
                #cax = plt.imshow(A, cmap='magma', interpolation='nearest')
                #fig.colorbar(cax)
                fig, ax = plt.subplots(figsize=(6,12))
                ax = sns.heatmap(A[0:249,0:23], cmap = 'viridis')
                #fig.show()
                fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/wpli/', e, '_', s, '_', b, '_mm_connection_matrix.png']), format = 'png')
                perf = A[-1,:]
                perf = np.reshape(perf, (1,perf.size))
                fig, ax = plt.subplots(figsize=(6,2))
                cmap = sns.cubehelix_palette(light=1, as_cmap=True)
                ax = sns.heatmap(perf, cmap = cmap)        
                fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/wpli/', e, '_', s, '_', b, '_perf.pdf']), format = 'pdf')



#%% Plot Dynamics FC

session = ['Sesion1', 'Session2', 'Session3', 'Session4']
tests = ['test06']
nNode = 102

for e in sens:
    for s in subjs:
        for b in bands: 
            for sess in session:
                for t in tests:
                    # Load
                    my_file = Path("".join(['/Users/stiso/Documents/MATLAB/NetBCI/', sess, '/', t, '/', s, '/FCmatrices/NMF_', b, '_', e, '_wpli.mat']))
                    if my_file.is_file():
                
                        # now plot subgraphs
                        dfc = sio.loadmat("".join(['/Users/stiso/Documents/MATLAB/NetBCI/', sess, '/', t, '/', s, '/FCmatrices/NMF_', b, '_', e, '_wpli.mat']))
                        A_vect = dfc.get('wpli_grad')
                        for trial in range(A_vect.shape[1]):
                            curr = A_vect[:,trial]
                            A = np.zeros([nNode,nNode])
                            cnt = -1
                            for col in range(nNode):
                                for row in range(col):
                                    if (row != col):
                                        cnt = cnt + 1;
                                        A[row,col] = curr[cnt]
                            A = A + np.transpose(A)
                
                            fig = plt.figure()
                            cax = plt.imshow(A, cmap='viridis', interpolation='nearest')
                            fig.colorbar(cax)
                            fig.show()
                            fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/wpli/', e, '_', s, '_', b, '_', sess, '_', t, '_', str(trial), '_dfc.pdf']), format = 'pdf')
                
                
                
                