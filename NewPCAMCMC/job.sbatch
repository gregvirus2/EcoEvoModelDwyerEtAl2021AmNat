#!/bin/bash

#SBATCH --job-name=GA8O42ZC
#SBATCH --output=GA8O42ZCOut.txt
#SBATCH --error=GA8O42ZCOuterror.txt 
#SBATCH --time=36:00:00
#SBATCH --array=1
#SBATCH --ntasks=1

module load openmpi


		

mpiexec ./NewPCAMCMC


