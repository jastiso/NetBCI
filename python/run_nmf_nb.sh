#!/bin/sh

subjects=( '001' '002' '003' '004' '005' '006' '007' '008' '009' '010' '011' '012' '013' '014' '015' '016' '017' '018' '019' '020')
bands=( 'alpha' 'beta' 'low_gamma' 'gamma')
sens=( 'grad')
for i in "${subjects[@]}"
do
for j in "${bands[@]}"
do
for k in "${sens[@]}"
do	
	qsub -o '/home/jstiso/output_files/' -e '/home/jstiso/error_files/' ./qsub_nmf_nb.sh $i $j $k
done
done
done
