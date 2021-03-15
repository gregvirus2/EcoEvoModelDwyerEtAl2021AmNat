#!/bin/bash

#SBATCH --job-name=PA_N1e5CNE8O42ZC
#SBATCH --output=t_WN1e5CNE8O42ZC.txt
#SBATCH --error=e_WN1e5CNE8O42ZCOuterror.txt 
#SBATCH --time=36:00:00
#SBATCH --array=1
#SBATCH --ntasks=1

module load openmpi

	

mpiexec ./DLPCAMCMC


