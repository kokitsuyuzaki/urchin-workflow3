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

if [ $1 = "sp5_0_GCF_genomic.fa" ]; then
    wget -P data/echinobase/ http://ftp.echinobase.org/pub/Genomics/Spur5.0/sp5_0_GCF_genomic.fa.gz --no-check-certificate
    cd data/echinobase
    gunzip sp5_0_GCF_genomic.fa.gz
else
    wget -P data/echinobase/ http://ftp.echinobase.org/pub/Genomics/Spur5.0/$1 --no-check-certificate
fi
