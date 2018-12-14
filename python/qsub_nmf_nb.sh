#!/bin/bash

subj=$1
band=$2
sens=$3
cd /data/jag/bassett-lab/jstiso/Python/NetBCI/

python ./NMF_pipe_gc_nb.py $subj $band $sens
python ./NMF_pipe_pr_nb.py $subj $band $sens
