#!/bin/sh

subjects=( '002')
bands=( 'theta' 'alpha' 'beta' 'gamma')
sens=( 'mag')
for i in "${subjects[@]}"
do
for j in "${bands[@]}"
do
for k in "${sens[@]}"
do	
	./qsub_bnv.sh $i $j $k
done
done
done
