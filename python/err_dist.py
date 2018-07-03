#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 14:17:19 2018

@author: stiso
"""

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
sens = ['mag', 'grad']


for e in sens:
    err = []
    err_bl = []
    for s in subjs:
        for b in bands:
            # Load
                my_file = Path("".join([top_dir, s, '/', e, '/gc_', b ,'_err.npy']))
                if my_file.is_file():
                    er = np.load("".join([top_dir, s, '/', e, '/gc_', b ,'_err.npy']))
                    err.append(er[er.size-1])

                               
        
#Now for Baseline
#    for s in subjs:
#        for b in bands:
#                my_file = Path("".join([top_dir, s, '/', e, '/ac_', b ,'_err_bl.npy']))
#                if my_file.is_file():
                    # Load
#                    er = np.load("".join([top_dir, s, '/', e, '/ac_', b ,'_err_bl.npy']))
#                    err_bl.append(er[er.size-1])
                
                
                # Plot
                
    fig = plt.figure()

    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(err, num_bins, normed=1, facecolor='orange', alpha=0.5)
    plt.xlabel('Error')
    plt.ylabel('Frequency')
    fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/gc/', e, '_error_hist.png']))


    # Plot baseline
                
    #fig = plt.figure()

    #num_bins = 50
    # the histogram of the data
    #n, bins, patches = plt.hist(err_bl, num_bins, normed=1, facecolor='grey', alpha=0.5)
    #plt.xlabel('Error')
    #plt.ylabel('Frequency')
    #fig.savefig("".join(['/Users/stiso/Documents/Python/NetBCI/NMF/ac/', e, '_bl_error_hist.png']))


