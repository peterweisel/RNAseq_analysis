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

module load miniconda

### assign directories
DATA=/projects/nomalab/shared/pweisel_project/Peter_data/Tbp1
RSEM_STAR=/projects/nomalab/pweisel/RSEM_STAR/pombe

### declare file names (files should be in a gzip fastq format)
declare -a lst=("TTAGGC_S7") # "TGACCA_S8" "ACAGTG_S9" "GCCAAT_S10" "ACTTGA_S11" "GATCAG_S12"

for NAME in ${lst[@]};
do
    # calculate expression
    conda activate rsem
    rsem-calculate-expression -p 12 --bam ${DATA}/${NAME}/${NAME}_Aligned.toTranscriptome.out.bam ${RSEM_STAR}/pombe ${NAME}
done