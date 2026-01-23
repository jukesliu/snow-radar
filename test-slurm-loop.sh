#!/bin/bash

for wtmm_scale in {800..1800..200}
do 
	echo $wtmm_scale
	for size_thresh in {50..200..50}
	do
		echo $size_thresh
		sbatch test-slurm.sh 500 $wtmm_scale $size_thresh 1
	done
done

