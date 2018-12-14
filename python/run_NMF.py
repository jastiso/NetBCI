#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 11:42:20 2018

@author: stiso
"""

import importlib.util
spec = importlib.util.spec_from_file_location('NMF_pipe_gc', '/data/jag/bassett-lab/jstiso/Python/NetBCI/NMF_pipe_gc.py')
NMF_pipe_gc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(NMF_pipe_gc)

spec = importlib.util.spec_from_file_location('NMF_pipe_ts', '/data/jag/bassett-lab/jstiso/Python/NetBCI/NMF_pipe_ts.py')
NMF_pipe_ts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(NMF_pipe_ts)


subjs = ['001', '002','003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019', '020']
bands = ['alpha', 'beta', 'low_gamma', 'gamma']
sens = ['grad', 'mag']

#%%

for s in subjs:
    for b in bands:
        for e in sens:
           NMF_pipe_gc.pipe(s,b,e)
           NMF_pipe_ts.pipe(s,b,e)
