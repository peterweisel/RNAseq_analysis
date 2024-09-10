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
# Description: This script uses UCSC utilities to calculate ChIP-seq binding frequencies for each gene in the specified expression groups. 
#              It is designed to run on the University of Oregon Talapas HPC cluster using SLURM job scheduler.
# Date: 01/20/2024

### load modules
module load ucsc-userapps/20191024

### assign directories
DIR_DATA=/projects/nomalab/pweisel/ChIPseq_analysis/data
DIR_REF=/projects/nomalab/pweisel/ChIPseq_analysis/genes/expression_groups
LOG_FILE=${DIR_DATA}/ave_score/log.txt 

### declare group and sample names
declare -a lst=("GROUP_1" "GROUP_2" "GROUP_3" "GROUP_4" "GROUP_5" "GROUP_6" "GROUP_7" "GROUP_8" "GROUP_9" "GROUP_10")
declare -a samples=("Cut14Pk")

for GRP in "${lst[@]}"; do
    for BW_NAME in "${samples[@]}"; do
        DIR_OUT="${DIR_DATA}/ave_score/${BW_NAME}/${GRP}"
        mkdir -p "${DIR_OUT}"

        FILE_LIST="${DIR_REF}/grouped_genes/${GRP}_genes.txt"
        GENE_LIST=$(cat "${FILE_LIST}")

        for GENE_NAME in ${GENE_LIST}; do 
            bigWigAverageOverBed -minMax ${DIR_DATA}/bigwig/${BW_NAME}.bw ${DIR_REF}/${GRP}/${GENE_NAME}.bed ${DIR_OUT}/${GENE_NAME}.txt
        done
    done
done

