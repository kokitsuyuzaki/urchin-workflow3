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

outdir=`echo $2 | sed -e 's|/plot/Allstates.png||'`

mkdir -p tmp2
cd tmp2
cp /landscaper .
cp /Snakefile .
cp -rf /src .

./landscaper -i ../$1 -o ../$outdir --cores=4 --memgb=10
