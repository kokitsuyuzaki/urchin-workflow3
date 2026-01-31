#!/bin/bash
#$ -l nc=4
#$ -p -50
#$ -r yes
#$ -q node.q

#SBATCH -n 4
#SBATCH --nice=50
#SBATCH --requeue
#SBATCH -p node03-06
#SBATCH --mem=32G
SLURM_RESTART_COUNT=2

BAMDIR="output/hpbase/$1/outs"
INFILE="${BAMDIR}/possorted_genome_bam.bam"
OUTFILE="${BAMDIR}/cellsorted_possorted_genome_bam.bam"  # この名前が重要！

# 壊れてる/未完なら作り直す
if [ -e "$OUTFILE" ]; then
  samtools quickcheck -q "$OUTFILE"
  if [ $? -ne 0 ]; then
    echo "[WARN] $OUTFILE is broken. Rebuilding..."
    rm -f "$OUTFILE"
  fi
fi

# 無いなら作る
if [ ! -e "$OUTFILE" ]; then
  samtools sort -t CB -O BAM -@ 4 -m 4G -o "$OUTFILE" "$INFILE"
fi

# ソート済みファイルが存在すればvelocytoはそれを使う
velocyto run10x "output/hpbase/$1" "data/hpbase/HpulGenome_v1_geneid.gtf"