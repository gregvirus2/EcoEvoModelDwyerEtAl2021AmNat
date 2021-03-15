#!/bin/sh
#PBS -N W_GA8O15ZC4
#PBS -S /bin/sh
#PBS -d .
#PBS -V


for (( i=0; i<$1;i++)); do
	echo `sbatch job.sbatch`	
done

