#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 11:42:20 2018

@author: stiso
"""

from NMF_pipe import pipe
from NMF_pipe_bl import pipe_bl

subjs = ['002','003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016']
bands = ['alpha', 'beta', 'low_gamma', 'gamma', 'all']
sens = ['grad', 'mag']

#%%
for s in subjs:
    for b in bands:
        for e in sens:
            pipe(s,b,e)
