#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 10:17:05 2018

@author: stiso
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 15:53:15 2018

@author: stiso
"""

import numpy as np
from pathlib import Path
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import plotly.plotly as py
from mpl_toolkits.axes_grid1 import make_axes_locatable

top_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/gc/'
subjs = ['001',  '002','003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018','019', '020']
bands = ['alpha', 'beta', 'low_gamma', 'gamma']
sens = ['grad']
nNode = 102

#%%
for e in sens:
    for s in subjs:
        for b in bands: 
            # Load
            my_file = Path("".join([top_dir, s, '/', e, '/nb_gc_', b ,'_coeff.npy']))
            if my_file.is_file():
                coeff = np.array(np.load("".join([top_dir, s, '/', e, '/nb_gc_', b ,'_coeff.npy'])))
                fig = plt.figure()
                for i in range(coeff.shape[0]):
                    plt.plot(coeff[i])
                fig.show()
                fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/gc/', e, '_', s, '_', b, '_coeffs_nb.png']))
                
                # now plot subgraphs
                subset = np.array(np.load("".join([top_dir, s, '/', e, '/nb_gc_', b ,'_subset.npy'])))
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
                    cax = plt.imshow(A, cmap='magma', interpolation='nearest')
                    fig.colorbar(cax)
                    fig.show()
                    if (subset[i,subset.shape[1]-1] == high):
                        fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/gc/', e, '_', s, '_', b, '_subset_high_nb.png']))
                    elif (subset[i,subset.shape[1]-1] == low):
                        fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/gc/', e, '_', s, '_', b, '_subset_low_nb.png']))
                    else:
                         fig.savefig(''.join(['/Users/stiso/Documents/Python/NetBCI/images/gc/', e, '_', s, '_', b, '_subset_', str(i), '_nb.png']))
                       
