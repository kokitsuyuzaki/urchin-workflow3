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

# Echinobase
cd data/echinobase
rm -rf sp5_0
cellranger mkref --genome=sp5_0 --fasta=../../$2 --genes=../../$1
