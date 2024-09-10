# !/usr/bin/Rscript
# Author: Peter Joseph Weisel
# Affiliation: Noma Laboratory at the University of Oregon
# Description: This R script processes ChIP-seq average scores, computes summary statistics for each gene based on expression level,
#              and generates a line plot depicting the average scores for different expression groups.
# Date: 01/20/2024

### load libraries
suppressWarnings({
  library(dplyr)
  library(data.table)
  library(pbapply)
  library(ggplot2)
  library(ggsci)
  library(cowplot)
  library(RColorBrewer)
  library(reshape2)
  library(reshape)
  library(tidyr)
  library(hrbrthemes)
  library(ggh4x)
  library(glue)
  library(scales)
})


options(scipen = 4)
R_table_all <- c()

### establish sample and expression group vectors
SAMPLES <- c("Cut14Pk")
GRP_LIST <- c("GROUP_1", "GROUP_2", "GROUP_3", "GROUP_4", "GROUP_5", "GROUP_6", "GROUP_7", "GROUP_8", "GROUP_9", "GROUP_10")


for (GRP in GRP_LIST){
  
  for (SCORE_NAME in SAMPLES){
    # establish directories
    DATA_DIR <- "/Users/peterweisel/Desktop/noma/ChIPseq/ave_score/"
    SCORE_DIR <- glue("/Users/peterweisel/Desktop/noma/ChIPseq/ave_score/{SCORE_NAME}/{GRP}/")
    FILE_gene <- fread(glue("/Users/peterweisel/Desktop/noma/genome/genes/expression_groups/grouped_genes/{GRP}_genes.txt"), head = FALSE, colClasses = 'character') 
    GENE_LIST <- as.list(t(FILE_gene))
    
    C_score <- c()
    for( SCORE_TXT in GENE_LIST){
      FILE_in <- paste0(SCORE_DIR, SCORE_TXT, ".txt")
      D_score <- fread(FILE_in)
      colnames(D_score)[1] <- "location"
      colnames(D_score)[5] <- "score"   # mean0
      D_score <- D_score %>%
        mutate( category = rep(SCORE_TXT, nrow(D_score))) %>%
        select(location, score, category)
      
      C_score <- rbind.data.frame(C_score, D_score)
    }
    
    R_table <- cast(C_score, category~location, value = "score")
    R_table[is.na(R_table)] <- 0
    R_table2 <- R_table %>%
      summarise_if(is.numeric, mean)
    R_table3 <- melt(R_table2)
    R_table3$score_name <- rep(SCORE_NAME, nrow(R_table3))
    R_table_all <- rbind.data.frame(R_table_all, R_table3)
    
  }
}

R_table_all_new_1 <- R_table_all
R_table_all_new_1$score_name <- factor(R_table_all_new_1$score_name,                 # Relevel group factor
                                       c("Cut14Pk"))
R_table_all_new_1 <- tibble::rowid_to_column(R_table_all_new_1, "ID")

sample_DF <- R_table_all_new_1
sample_DF$group <- as.factor(ifelse(sample_DF$ID<201, 'GROUP_1',
                                    ifelse(sample_DF$ID<401, 'GROUP_2',
                                           ifelse(sample_DF$ID<601, 'GROUP_3',
                                                  ifelse(sample_DF$ID<801, 'GROUP_4',
                                                         ifelse(sample_DF$ID<1001, 'GROUP_5',
                                                                ifelse(sample_DF$ID<1201, 'GROUP_6',
                                                                       ifelse(sample_DF$ID<1401, 'GROUP_7',
                                                                              ifelse(sample_DF$ID<1601, 'GROUP_8',
                                                                                     ifelse(sample_DF$ID<1801, 'GROUP_9',
                                                                                            ifelse(sample_DF$ID<2001, 'GROUP_10')))))))))))

sample_DF$group <- factor(sample_DF$group, levels = c("GROUP_1", "GROUP_2", "GROUP_3", "GROUP_4", "GROUP_5", "GROUP_6", "GROUP_7", "GROUP_8", "GROUP_9", "GROUP_10"))

p1 <- ggplot(sample_DF, aes(x= as.numeric(variable), y= value, group=group, color=group)) +
  geom_line(linewidth=0.75, key_glyph = "rect") +
  scale_color_manual("Expression\nLevel", labels = c("High", "", "", "", "", "", "", "", "", "Low"), values=c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')) +
  scale_x_continuous(minor_breaks = seq(0, 200, by = 5), breaks = c(0, 25, 50, 75, 100, 105, 110, 115, 120, 125, 150, 175, 200), labels = c("-2kb", "-1.5kb", "-1kb", "-0.5kb", "0%", "20%", "40%", "60%", "80%", "100%", "+1kb", "+1.5kb", "+2kb"), guide = "axis_minor") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(
    
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    legend.text=element_text(size= 10, angle = 0, hjust = 0),
    legend.spacing.y = unit(0.25, 'cm'),
    legend.position= "right",
    legend.margin=margin(t=-10),
    legend.title=element_text(size=10),
    plot.title = element_text(hjust = 0, size=10, vjust = 0) #, ,
  ) +
  labs(x="Distance From TSS", y="Average Cut14Pk Score", title="Cut14Pk Score for Pol II genes")

p1
