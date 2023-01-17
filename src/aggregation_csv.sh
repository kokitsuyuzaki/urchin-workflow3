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

echo "sample_id,molecule_h5" > $1
echo "cont-24h,cont-24h/outs/molecule_info.h5" >> $1
echo "cont-48h,cont-48h/outs/molecule_info.h5" >> $1
echo "cont-72h,cont-72h/outs/molecule_info.h5" >> $1
echo "cont-96h,cont-96h/outs/molecule_info.h5" >> $1
echo "DAPT-24h,DAPT-24h/outs/molecule_info.h5" >> $1
echo "DAPT-48h,DAPT-48h/outs/molecule_info.h5" >> $1
echo "DAPT-72h,DAPT-72h/outs/molecule_info.h5" >> $1
echo "DAPT-96h,DAPT-96h/outs/molecule_info.h5" >> $1
