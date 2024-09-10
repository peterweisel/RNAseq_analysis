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
# Description: This script aligns Paired-end ChIP-seq FASTQ files to the Bowtie2 reference genome of Schizosaccharomyces pombe.
#              It is designed to run on the University of Oregon Talapas HPC cluster using SLURM job scheduler.
# Date: 01/20/2024

### load modules
module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load easybuild  ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load eb-hide/1  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load eb-hide/1  ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load SAMtools/1.8
module load miniconda
module load bowtie2/2.4.4

### assign directories 
DIR_WORK=/projects/nomalab/pweisel/ChIPseq_analysis
DIR_DATA=${DIR_WORK}/Data_raw

### declare sample names
declare -a lst=("Cut14Pk")
for NAME in ${lst[@]};
do

    # output files will be placed in the following directory
    DIR_OUT=${DIR_WORK}/alignment_files/${NAME}
    mkdir -p ${DIR_OUT}
    cd ${DIR_OUT}

    # pired-end ChIP-seq FASTQ files
    FASTQ_FILE_1=${NAME}_pass1.fastq
    FASTQ_FILE_2=${NAME}_pass2.fastq

    # display compressed FASTQ files
    zcat < ${DIR_DATA}/${FASTQ_FILE_1}.gz > ${FASTQ_FILE_1}
    zcat < ${DIR_DATA}/${FASTQ_FILE_2}.gz > ${FASTQ_FILE_2}
    
    # Bowtie2 alignment
    bowtie2 -p 4 --no-unal -x ${DIR_WORK}/genome/data/pombe/2018/Bowtie2/pombe -1 ${FASTQ_FILE_1} -2 ${FASTQ_FILE_2} -S ${NAME}.sam 2> ${NAME}.log

    # sort and index BAM file
    samtools view -f 0x2 -h -b ${NAME}.sam > ${NAME}.bam
    samtools sort ${NAME}.bam > ${NAME}_sorted.bam
    samtools index ${NAME}_sorted.bam

    # perform normalization and generate genome coverage
    conda activate deeptools-3.5.1  
    bamCoverage --ignoreForNormalization III --normalizeUsing RPGC --effectiveGenomeSize 10118937 -b ${NAME}_sorted.bam -o ${NAME}_normalized.bw

done
