subj=$1
band=$2
sens=$3
cd /data/jag/bassett-lab/jstiso/Python/NetBCI/

python ./NMF_pipe_wpli.py $subj $band $sens

