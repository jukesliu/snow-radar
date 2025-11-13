#!/bin/bash

for wtmm_scale in {600..1000..200}
do 
	echo $wtmm_scale
	for size_thresh in {50..150..50}
	do
		echo $size_thresh
		sbatch test-slurm.sh 800 $wtmm_scale $size_thresh 1
	done
done

