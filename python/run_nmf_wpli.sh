#!/bin/sh

#subjects=( '001' '002' '003' '004' '005' '006' '007' '008' '009' '010' '011' '012' '013' '014' '015' '016' '017' '018' '019' '020')
subjects=( '017' '010')
bands=( 'beta' 'low_gamma')
sens=( 'grad')
for i in "${subjects[@]}"
do
for j in "${bands[@]}"
do
for k in "${sens[@]}"
do	
	qsub -o '/home/jstiso/output_files/' -e '/home/jstiso/error_files/' ./qsub_nmf_wpli.sh $i $j $k
	#qsub -o '/home/jstiso/output_files/' -e '/home/jstiso/error_files/' ./qsub_nmf_wpli_pr.sh $i $j $k
	#qsub -o '/home/jstiso/output_files/' -e '/home/jstiso/error_files/' ./qsub_nmf_wpli_ind.sh $i $j $k
done
done
done
