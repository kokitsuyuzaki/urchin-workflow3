# #!/bin/bash
# #$ -l nc=4
# #$ -p -50
# #$ -r yes
# #$ -q node.q

# #SBATCH -n 4
# #SBATCH --nice=50
# #SBATCH --requeue
# #SBATCH -p node03-06
# SLURM_RESTART_COUNT=2

# outdir=`echo $2 | sed -e 's|/BIN_DATA.tsv||'`

# mkdir -p tmp
# cd tmp
# cp /sbmfcv .
# cp /Snakefile .
# cp -rf /src .

# ./sbmfcv -i ../$1 -o ../$outdir \
# --cores=8 --rank_min=5 --rank_max=10 \
# --lambda_min=-3 --lambda_max=3 --trials=50 \
# --n_iter_max=100 --ratio=20 --memgb=50

# rm -rf tmp

# 単に値を-1,1に変換
awk 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $i=($i>0)?1:($i<0)?-1:0; print}' $1 > $2