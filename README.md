# Roles-of-TBP-N-Terminus-in-3D-Genome-Organization-and-Gene-Regulation

The general transcription factor, TATA box-binding protein (TBP), plays a crucial role in gene transcription by converting genetic information from DNA to RNA. Although this factor is one of the most important transcription factors involved in every type of gene expression and is evolutionarily conserved in eukaryotes, the role of TBP in three-dimensional (3D) genome organization remains unclear. To this end, our recent research has highlighted the multifaceted functions of the TBP N-terminus in transcriptional regulation and 3D genome organization. Here, we show that the N-terminus deletion of TBP (TBP∆N) impairs cell growth at several culturing conditions. Also, our genomics analysis using next-generation sequencing technology reveals that a specific 3D genome structure(gene-sized small chromatin domains), is disrupted in the TBP∆N mutants. Furthermore, we show that the TBP N∆ mutants display a global suppression of gene expression, although there were no significant changes in protein binding of wild-type and mutant TBP proteins across the fission yeast genome. Considering these findings, we propose that the TBP N-terminus is an essential factor in producing proper stress responses in fission yeast cells, which is critical for maintaining faithful segregation of transcribed genes, in addition to their expression and localization across the genome.

## Overview

This repository contains scripts for the analysis of RNA-seq data for the Tbp1 protein in fission yeast. The analysis includes STAR alignment, RSEM quantification, and differential gene expression using edgeR.

## Scripts

1. edgeR_practice.R: Differential gene expression analysis using edgeR with custom dispersion value.

2. rsem.sh: Estimate gene and isoform expression from RNA-Seq data.

3. star.sh: Alignment of raw FASTQ files to S. pombe reference genome using STAR and samtools.
