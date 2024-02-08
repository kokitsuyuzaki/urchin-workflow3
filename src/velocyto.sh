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

SAMPLE="output/hpbase/"$1
GTF="data/hpbase/HpulGenome_v1_geneid.gtf"

samtools --version
INFILE="output/hpbase/"$1"/outs/possorted_genome_bam.bam"
OUTFILE="output/hpbase/"$1"/outs/cellsorted_possorted_genome_bam.bam"
if [ ! -e $OUTFILE ]; then
	samtools sort -t CB -O BAM -o $OUTFILE $INFILE
fi

velocyto run10x $SAMPLE $GTF
