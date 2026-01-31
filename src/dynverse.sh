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

export APPTAINER_NO_HOME=1
export SINGULARITY_NO_HOME=1
export APPTAINER_CONTAINALL=1
export SINGULARITY_CONTAINALL=1

Rscript src/dynverse.R $@
