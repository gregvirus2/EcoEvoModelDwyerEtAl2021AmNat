#!/bin/sh
#PBS -N PA_N1e5CNE8O42ZC
#PBS -S /bin/sh
#PBS -d .
#PBS -V


for (( i=0; i<$1;i++)); do
	echo `sbatch job.sbatch`	
done

