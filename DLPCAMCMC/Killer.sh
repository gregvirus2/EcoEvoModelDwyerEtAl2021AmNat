#!/bin/bash
#$ -S /bin/bash
#$ -cwd


for (( i=$1; i<$1+$2;i++)); do
        	echo `qdel $i`
done