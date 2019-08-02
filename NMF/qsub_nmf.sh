#!/bin/bash

subj=$1
band=$2
sens=$3
cd /data/jag/bassett-lab/jstiso/Python/NetBCI/

python ./NMF_pipe_gc.py $subj $band $sens
python ./NMF_pipe_pr.py $subj $band $sens
python ./NMF_pipe_ind.py $subj $band $sens
