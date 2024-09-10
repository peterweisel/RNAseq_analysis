#!/bin/bash

### align
WORK=/Users/nomalab/Projects/Peter
RSEM_STAR=/Users/nomalab/Projects/Peter/Genome/Saccharomyces_cerevisiae/R64-1-1/RSEM_star

# pombe = /Users/nomalab/Projects/Peter/Genome/pombe/RSEM_STAR
# cerevisiae = /Users/nomalab/Projects/Peter/Genome/Saccharomyces_cerevisiae/R64-1-1/RSEM_star


declare -a lst=("ACAGTG_S9_001" "ACTTGA_S11_001" "GATCAG_S12_001" "TGACCA_S8_001" "TTAGGC_S7_001" "GCCAAT_S10_001")
for NAME in ${lst[@]};
do
    OUTPUT=/Users/nomalab/Projects/Peter/Data/RNA-seq/tbp1_cerevisiae/${NAME}
    mkdir -p ${OUTPUT}
    cd ${OUTPUT}

    gzip -d ${WORK}/Data_raw/${NAME}_pass_1.fastq.gz
    gzip -d ${WORK}/Data_raw/${NAME}_pass_2.fastq.gz

    STAR --genomeDir ${RSEM_STAR} \
        --runThreadN 8 \
        --readFilesIn ${WORK}/Data_raw/${NAME}_pass_1.fastq,${WORK}/Data_raw/${NAME}_pass_2.fastq \
        --outSAMattributes MD \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix "${OUTPUT}/${NAME}_" \
        -outSAMtype BAM SortedByCoordinate

    gzip ${WORK}/Data_raw/${NAME}_pass_1.fastq
    gzip ${WORK}/Data_raw/${NAME}_pass_2.fastq

    samtools view -b ${NAME}_Aligned.out.sam > ${NAME}.bam

    rm ${NAME}_Aligned.out.sam

    # Sort by chromosome number to make a bigwig file.
    # samtools sort ${NAME}.bam > ${NAME}_sorted.bam

    # Make an index file of the bam file.
    samtools index ${NAME}.bam

    ### make bam index file

    ### bamCoverage
    source /Users/nomalab/opt/anaconda3/etc/profile.d/conda.sh
    conda activate deeptools
    bamCoverage --numberOfProcessors max --binSize 10 -b ${NAME}.bam -o ${NAME}_Aligned.sortedByCoord.out.bw 2> ${NAME}_Aligned.sortedByCoord.out.bw.log

done
