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

if [ $1 = "hpbase" ]; then
	SAMPLE="output/hpbase/"$2
	GTF="data/hpbase/HpulGenome_v1_geneid.gtf"
else
	SAMPLE="output/echinobase/"$2
	GTF="data/echinobase/sp5_0_GCF_geneid.gtf"
fi

samtools --version
INFILE="output/"$1"/"$2"/outs/possorted_genome_bam.bam"
OUTFILE="output/"$1"/"$2"/outs/cellsorted_possorted_genome_bam.bam"
if [ ! -e $OUTFILE ]; then
	samtools sort -t CB -O BAM -o $OUTFILE $INFILE
fi

velocyto run10x $SAMPLE $GTF
