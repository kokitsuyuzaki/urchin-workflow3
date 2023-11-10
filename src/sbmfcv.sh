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

outdir=`echo $2 | sed -e 's|/BIN_DATA.tsv||'`

mkdir -p tmp
cd tmp
cp /sbmfcv .
cp /Snakefile .
cp -rf /src .

./sbmfcv -i ../$1 -o ../$outdir \
--cores=8 --rank_min=2 --rank_max=15 \
--lambda_min=-3 --lambda_max=3 --trials=5 \
--n_iter_max=50 --ratio=20 --memgb=50

rm -rf tmp
