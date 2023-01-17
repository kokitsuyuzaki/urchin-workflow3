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
cd data/hpbase
rm -rf HpulGenome_v1
cellranger mkref --genome=HpulGenome_v1 --fasta=../../$2 --genes=../../$1
