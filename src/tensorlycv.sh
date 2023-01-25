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

outdir=`echo $2 | sed -e 's|/plot/bestrank_besttrial_barplot_FINISH||'`

mkdir -p tmp
cd tmp
cp /tensorlycv .
cp /Snakefile .
cp -rf /src .

./tensorlycv -i ../$1 -o ../$outdir \
--cores=4 --rank=8 --trials=50 --iters=1000 \
--ratio=1 --memgb=10

cd ..
rm -rf tmp