#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=RNAseq_tbp1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pweisel@uoregon.edu
#SBATCH --error=job.%J.err      
#SBATCH --time=0-23:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=nomalab

# Author: Peter Joseph Weisel
# Affiliation: Noma Laboratory at the University of Oregon
# Description: This script uses MACS2 to perform peak calling on ChIP-seq treatment and control samples.
#              It is designed to run on the University of Oregon Talapas HPC cluster using SLURM job scheduler.
# Date: 01/20/2024

### load modules
module load racs-eb/1
module load MACS2/2.1.1.20160309-intel-2017b-Python-2.7.14

### directory with alignment files
DIR_DATA=/projects/nomalab/pweisel/ChIPseq_analysis/alignment_files

declare -a lst=("Cut14Pk")
for NAME_TARGET in ${lst[@]};
do

    # control sample
    NAME_CONTROL=${NAME_TARGET}_control

    # peak calling files will be placed in the following directory
    DIR_OUT=${DIR_DATA}/peak_calling/${NAME_TARGET}
    mkdir -p ${DIR_OUT}

    # chromosome length of S.pombe
    CHROM_LENGTH=12571820

    # perform peak calling
    macs2 callpeak -t ${DIR_DATA}/${NAME_TARGET}/${NAME_TARGET}_sorted.bam -c ${DIR_DATA}/${NAME_CONTROL}/${NAME_CONTROL}_sorted.bam -g ${CHROM_LENGTH} --bdg -f BAM -n ${NAME_TARGET} --nomodel --extsize 147 --outdir ${DIR_OUT} --keep-dup all 2> ${NAME_TARGET}_macs2_callpeak.log

done