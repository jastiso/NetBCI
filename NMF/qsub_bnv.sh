#!/bin/bash

subj=$1
band=$2
nGraph=$4
sens=$3
cd /data/jag/bassett-lab/jstiso/Python/NetBCI/

matlab -r "plot_bnv('$subj', '$band', 6, '$sens'); exit"
