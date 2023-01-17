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

rm -rf output/hpbase/$1
mkdir output/hpbase/$1
cd output/hpbase

cellranger count \
    --id=$1 \
    --transcriptome=../../data/hpbase/HpulGenome_v1 \
    --fastqs=../../data/Azenta/00_Rawdata/$1 \
    --sample=$1 \
    --localcores=4 \
    --localmem=64
