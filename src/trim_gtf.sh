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
cat $1 | awk -F ' ' '{if($11 == "gene_id") print}' > $3

# Echinobase
cat $2 | awk -F ' ' '{if($11 == "gene_id") print}' > $4
