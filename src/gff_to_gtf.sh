#!/bin/bash
#$ -l nc=4
#$ -p -50
#$ -r yes
#$ -q node.q

#SBATCH -n 4
#SBATCH --nice=50
#SBATCH --requeue
#SBATCH -p node03-06
SLURM_RESTART_COUNT=2

# Hpbase
gffread $1 -g $3 -E -T -o $5

# Echinobase
gffread $2 -g $4 -E -T -o $6
