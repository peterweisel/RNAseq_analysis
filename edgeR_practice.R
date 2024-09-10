# RNA-seq Differential Expression Analysis
# Peter Weisel
# 02-13-2024

library(edgeR)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(ggrepel)

# load RSEM files (contains expected counts)
counts_input <- read.delim("/Users/peterweisel/Desktop/noma/RNAseq/new_pombe/RSEM_STAR/Tbp1/input_counts.tsv", header = TRUE)
counts_biotin <- read.delim("/Users/peterweisel/Desktop/noma/RNAseq/new_pombe/RSEM_STAR/Tbp1/biotin_counts.tsv", header = TRUE)

# check dimensionality and first few rows
dim(counts_input)
head(counts_input)
str(counts_input)
dim(counts_biotin)
head(counts_biotin)
str(counts_biotin)

###### BIOTIN
### make dgList
sample_info <- c("WT", "N25", "N50")
sample_info <- factor(sample_info, levels = c("WT", "N25", "N50"))
dgList <- DGEList(counts = counts_biotin, group = factor(sample_info))
keep <- filterByExpr(y = dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
dgList2 <- dgList

### take standard dev
data_frame <- data.frame(gene_id = dgList$genes, samples = dgList$counts)
data_frame$std_dev <- apply(data_frame[, 2:4], 1, sd)
ordered <- data_frame[order(data_frame$std_dev), ]
max_val <- length(ordered$std_dev) * 0.90
top_90 <- ordered[1:max_val, ]

### calculate dispersion value
housekeeping <- select(top_90, samples.WT:samples.N50)
dgList_filtered <- DGEList(counts = housekeeping)
dgList_filtered$samples$group <- 1
disp_filtered <- estimateCommonDisp(dgList_filtered)
dispersion_val_biotin <- disp_filtered$common.dispersion
dispersion_val_biotin

###### INPUT
### make dgList
sample_info <- c("WT", "N25", "N50")
sample_info <- factor(sample_info, levels = c("WT", "N25", "N50"))
dgList <- DGEList(counts = counts_input, group = factor(sample_info))
keep <- filterByExpr(y = dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
dgList2 <- dgList

### take standard dev
data_frame <- data.frame(gene_id = dgList$genes, samples = dgList$counts)
data_frame$std_dev <- apply(data_frame[, 2:4], 1, sd)
ordered <- data_frame[order(data_frame$std_dev), ]
max_val <- length(ordered$std_dev) * 0.90
top_90 <- ordered[1:max_val, ]

### calculate dispersion value
housekeeping <- select(top_90, samples.WT:samples.N50)
dgList_filtered <- DGEList(counts = housekeeping)
dgList_filtered$samples$group <- 1
disp_filtered <- estimateCommonDisp(dgList_filtered)
dispersion_val_input <- disp_filtered$common.dispersion
dispersion_val_input

### N50 Input
sample_info <- c("WT", "N25", "N50")
sample_info <- factor(sample_info, levels = c("WT", "N50", "N25"))
dgList <- DGEList(counts = counts_input, group = factor(sample_info))
keep <- filterByExpr(y = dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
dgList2 <- dgList

# spike-in normalization, extract logFC for the dispersion of 0.029
scaling_factors <- c(WT = 2.645712623, N50 = 1.636098886, N25 = 2.2516944)
dgList$samples$lib.size <- colSums(dgList$counts)
dgList$samples$norm.factors <- scaling_factors
dgList <- estimateDisp(dgList)
et <- exactTest(object = dgList, dispersion = 0.1^6)
top_degs = topTags(object = et, n = "Inf")
summary(decideTests(object = et, lfc = 1))

# remove logCPM values < 1
N50_input_df <- data.frame(top_degs)
N50_input_df <- N50_input_df[N50_input_df$logCPM >= 1, ]

# smear plot
ggplot(N50_input_df, aes(x = logCPM, y=logFC, fill=FDR < 0.05, color=FDR < 0.05)) + 
  geom_point(pch=21, alpha=0.75) + theme_light() + ylab('log2(Tbp1∆50 Input / WT)') + 
  scale_fill_manual(values = alpha(c("black", "red"))) +
  scale_color_manual(values = alpha(c("black", "red"))) +
  ylim(-12,12)


head(N50_input_df[order(N50_input_df$logFC), ])
head(N50_input_df[order(-N50_input_df$logFC), ])

# significant genes
up_N50_input <- N50_input_df$gene_id[N50_input_df$logFC > 0 & N50_input_df$FDR < 0.05]
down_N50_input <- N50_input_df$gene_id[N50_input_df$logFC < 0 & N50_input_df$FDR < 0.05]

top_up_N50_input <- N50_input_df[N50_input_df$gene_id %in% up_N50_input, ]
top_down_N50_input <- N50_input_df[N50_input_df$gene_id %in% down_N50_input, ]

up_N50_input_gene <- top_up_N50_input$gene_id
down_N50_input_gene <- top_down_N50_input$gene_id

# classify genes based on expression level
N50_input_df <- N50_input_df %>%  
  mutate(Expression = ifelse(N50_input_df$gene_id %in% top_up_N50_input$gene_id, "Upregulated", 
                             ifelse(N50_input_df$gene_id %in% top_down_N50_input$gene_id, "Downregulated", "Neither")))

# volcano plot showing differentially expressed genes
ggplot(data=N50_input_df, aes(x=logFC, y=-log10(FDR), label=gene_id, color=Expression)) +
  scale_color_manual(values = c("Upregulated" = "deepskyblue3", "Downregulated" = "red", "Neither" = "black")) +
  ggtitle("Differential Expression of Input Tbp1 WT vs Input Tbp1∆50") +
  ylab("-log10(FDR)") + xlab("log2(Tbp1∆50/Tbp1)") +
  geom_point() + 
  geom_text_repel(color = "black", size = 4, max.overlaps=1) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype=3) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype=3) +
  theme_classic() 


### N25 Input
sample_info <- c("WT", "N25", "N50")
sample_info <- factor(sample_info, levels = c("WT", "N25", "N50"))
dgList <- DGEList(counts = counts_input, group = factor(sample_info))
dgList
dgList$samples
head(dgList$counts)
head(dgList$genes)
keep <- filterByExpr(y = dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]

# spike-in normalization
scaling_factors <- c(WT = 2.645712623, N25 = 2.2516944, N50 = 1.636098886)
dgList$samples$lib.size <- colSums(dgList$counts)
dgList$samples$norm.factors <- scaling_factors

dgList <- estimateDisp(dgList)
et <- exactTest(object = dgList, dispersion = 0.1^6)
top_degs = topTags(object = et, n = "Inf")
head(top_degs)
summary(decideTests(object = et, lfc = 1))

# remove logCPM values < 1
N25_input_df <- data.frame(top_degs)
N25_input_df <- N25_input_df[N25_input_df$logCPM >= 1, ]

# smear plot
ggplot(N25_input_df, aes(x = logCPM, y=logFC, fill=FDR < 0.05, color=FDR < 0.05)) + 
  geom_point(pch=21, alpha=0.75) + theme_light() + ylab('log2(Tbp1∆25 Input / WT)') + 
  scale_fill_manual(values = alpha(c("black", "red"))) +
  scale_color_manual(values = alpha(c("black", "red"))) +
  ylim(-12,12)

head(N25_input_df[order(N25_input_df$logFC), ])
head(N25_input_df[order(-N25_input_df$logFC), ])

# significant genes
up_N25_input <- N25_input_df$gene_id[N25_input_df$logFC > 0 & N25_input_df$FDR < 0.05]
down_N25_input <- N25_input_df$gene_id[N25_input_df$logFC < 0 & N25_input_df$FDR < 0.05]

top_up_N25_input <- N25_input_df[N25_input_df$gene_id %in% up_N25_input, ]
top_down_N25_input <- N25_input_df[N25_input_df$gene_id %in% down_N25_input, ]

up_N25_input_gene <- top_up_N25_input$gene_id
down_N25_input_gene <- top_down_N25_input$gene_id

# classify genes based on expression level
N25_input_df <- N25_input_df %>%  
  mutate(Expression = ifelse(N25_input_df$gene_id %in% top_up_N25_input$gene_id, "Upregulated", 
                             ifelse(N25_input_df$gene_id %in% top_down_N25_input$gene_id, "Downregulated", "Neither")))

# volcano plot showing differentially expressed genes
ggplot(data=N25_input_df, aes(x=logFC, y=-log10(FDR), label=gene_id, color=Expression)) +
  scale_color_manual(values = c("Upregulated" = "deepskyblue3", "Downregulated" = "red", "Neither" = "black")) +
  ggtitle("Differential Expression of Input Tbp1 WT vs Input Tbp1∆25") +
  ylab("-log10(FDR)") + xlab("log2(Tbp1∆25/Tbp1)") +
  geom_point() + 
  geom_text_repel(color = "black", size = 4, max.overlaps=1) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype=3) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype=3) +
  theme_classic()

### N50 Biotin
sample_info <- c("WT", "N25", "N50")
sample_info <- factor(sample_info, levels = c("WT", "N50", "N25"))
dgList <- DGEList(counts = counts_biotin, group = factor(sample_info))
dgList
dgList$samples
head(dgList$counts)
head(dgList$genes)
keep <- filterByExpr(y = dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]

# spike-in normalization
scaling_factors <- c(WT = 1.673780233, N50 = 1.205850788, N25 = 1.672996169)
dgList$samples$lib.size <- colSums(dgList$counts)
dgList$samples$norm.factors <- scaling_factors

dgList <- estimateDisp(dgList)
et <- exactTest(object = dgList, dispersion = 0.0029)
top_degs = topTags(object = et, n = "Inf")
head(top_degs)
summary(decideTests(object = et, lfc = 1))

# remove logCPM values < 1
N50_biotin_df <- data.frame(top_degs)
N50_biotin_df <- N50_biotin_df[N50_biotin_df$logCPM >= 1, ]

# smear plot
ggplot(N50_biotin_df, aes(x = logCPM, y=logFC, fill=FDR < 0.05, color=FDR < 0.05)) + 
  geom_point(pch=21, alpha=0.5) + theme_light() + ylab('log2(Tbp1∆50 Biotin / WT)') + 
  scale_fill_manual(values = alpha(c("black", "darkred"))) +
  scale_color_manual(values = alpha(c("black", "darkred")))

head(N50_biotin_df[order(N50_biotin_df$logFC), ])
head(N50_biotin_df[order(-N50_biotin_df$logFC), ])

# significant genes
up_N50_biotin <- N50_biotin_df$gene_id[N50_biotin_df$logFC > 0 & N50_biotin_df$FDR < 0.05]
down_N50_biotin <- N50_biotin_df$gene_id[N50_biotin_df$logFC < 0 & N50_biotin_df$FDR < 0.05]

top_up_N50_biotin <- N50_biotin_df[N50_biotin_df$gene_id %in% up_N50_biotin, ]
top_down_N50_biotin <- N50_biotin_df[N50_biotin_df$gene_id %in% down_N50_biotin, ]

up_N50_biotin_gene <- top_up_N50_biotin$gene_id
down_N50_biotin_gene <- top_down_N50_biotin$gene_id

# classify genes based on expression level
N50_biotin_df <- N50_biotin_df %>%  
  mutate(Expression = ifelse(N50_biotin_df$gene_id %in% top_up_N50_biotin$gene_id, "Upregulated", 
                             ifelse(N50_biotin_df$gene_id %in% top_down_N50_biotin$gene_id, "Downregulated", "Neither")))

# volcano plot showing differentially expressed genes
ggplot(data=N50_biotin_df, aes(x=logFC, y=-log10(FDR), label=gene_id, color=Expression)) +
  scale_color_manual(values = c("Upregulated" = "deepskyblue3", "Downregulated" = "darkred", "Neither" = "black")) +
  ggtitle("Differential Expression of Biotin Tbp1 WT vs Biotin Tbp1∆50") +
  ylab("-log10(FDR)") + xlab("log2(Tbp1∆50/Tbp1)") +
  geom_point() + 
  geom_text_repel(color = "black", size = 4, max.overlaps=1) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype=3) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype=3) +
  theme_classic() 

### N25 Biotin
sample_info <- c("WT", "N25", "N50")
sample_info <- factor(sample_info, levels = c("WT", "N25", "N50"))
dgList <- DGEList(counts = counts_biotin, group = factor(sample_info))
dgList
dgList$samples
head(dgList$counts)
head(dgList$genes)
keep <- filterByExpr(y = dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]

# spike-in normalization
scaling_factors <- c(WT = 1.673780233, N25 = 1.672996169, N50 = 1.205850788)
dgList$samples$lib.size <- colSums(dgList$counts)
dgList$samples$norm.factors <- scaling_factors

dgList <- estimateDisp(dgList)
et <- exactTest(object = dgList, dispersion = 0.0029)
top_degs = topTags(object = et, n = "Inf")
head(top_degs)
summary(decideTests(object = et, lfc = 1))

# remove logCPM values < 1
N25_biotin_df <- data.frame(top_degs)
N25_biotin_df <- N25_biotin_df[N25_biotin_df$logCPM >= 1, ]

# smear plot
ggplot(N25_biotin_df, aes(x = logCPM, y=logFC, fill=FDR < 0.05, color=FDR < 0.05)) + 
  geom_point(pch=21, alpha=0.5) + theme_light() + ylab('log2(Tbp1∆25 Biotin / WT)') + 
  scale_fill_manual(values = alpha(c("black", "darkred"))) +
  scale_color_manual(values = alpha(c("black", "darkred")))

head(N25_biotin_df[order(N25_biotin_df$logFC), ])
head(N25_biotin_df[order(-N25_biotin_df$logFC), ])

# significant genes
up_N25_biotin <- N25_biotin_df$gene_id[N25_biotin_df$logFC > 0 & N25_biotin_df$FDR < 0.05]
down_N25_biotin <- N25_biotin_df$gene_id[N25_biotin_df$logFC < 0 & N25_biotin_df$FDR < 0.05]

top_up_N25_biotin <- N25_biotin_df[N25_biotin_df$gene_id %in% up_N25_biotin, ]
top_down_N25_biotin <- N25_biotin_df[N25_biotin_df$gene_id %in% down_N25_biotin, ]

up_N25_biotin_gene <- top_up_N25_biotin$gene_id
down_N25_biotin_gene <- top_down_N25_biotin$gene_id

# classify genes based on expression level
N25_biotin_df <- N25_biotin_df %>%  
  mutate(Expression = ifelse(N25_biotin_df$gene_id %in% top_up_N25_biotin$gene_id, "Upregulated", 
                             ifelse(N25_biotin_df$gene_id %in% top_down_N25_biotin$gene_id, "Downregulated", "Neither")))
N25_biotin_df2 <- N25_biotin_df %>% filter(gene_id != "SPRRNA.45")
# volcano plot showing differentially expressed genes
ggplot(data=N25_biotin_df2, aes(x=logFC, y=-log10(FDR), label=gene_id, color=Expression)) +
  scale_color_manual(values = c("Upregulated" = "deepskyblue3", "Downregulated" = "darkred", "Neither" = "black")) +
  ggtitle("Differential Expression of Biotin Tbp1 WT vs Biotin Tbp1∆25") +
  ylab("-log10(FDR)") + xlab("log2(Tbp1∆25/Tbp1)") +
  geom_point() + 
  geom_text_repel(color = "black", size = 4, max.overlaps=1) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype=3) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype=3) +
  theme_classic() 

dev.off(dev.list()["RStudioGD"])

### venn diagram to show overlapping genes

## downregulated genes
downregulated_list <- list(
  "Tbp1∆50 Input" = down_N50_input_gene,
  "Tbp1∆25 Input" = down_N25_input_gene,
  "Tbp1∆50 Biotin" = down_N50_biotin_gene,
  "Tbp1∆25 Biotin" = down_N25_biotin_gene
)

myCol <- brewer.pal(4, "Pastel1")

venn.plot <- venn.diagram(
  x = downregulated_list,
  category.names = names(downregulated_list),
  filename = NULL,
  output = TRUE,
  fill = myCol,
  fontfamily = "sans",
  cat.fontfamily = "sans",
)

venn.plot <- plot(venn.plot)
grid.draw(venn.plot)

# shared downregulated genes across all 4 samples
down_shared <- Reduce(intersect, list(down_N50_input_gene, down_N50_biotin_gene))
print(down_shared)

dev.off(dev.list()["RStudioGD"])

## upregulated genes
upregulated_list
upregulated_list <- list(
  "Tbp1∆50 Input" = up_N50_input_gene,
  "Tbp1∆25 Input" = up_N25_input_gene,
  "Tbp1∆50 Biotin" = up_N50_biotin_gene,
  "Tbp1∆25 Biotin" = up_N25_biotin_gene
)

myCol <- brewer.pal(4, "Pastel1")

venn.plot <- venn.diagram(
  x = upregulated_list,
  category.names = names(upregulated_list),
  filename = NULL,
  output = TRUE,
  fill = myCol,
  fontfamily = "sans",
  cat.fontfamily = "sans",
)

venn.plot <- plot(venn.plot)
grid.draw(venn.plot)

# shared upregulated genes across all 4 samples
up_shared <- Reduce(intersect, list(up_N50_input_gene, up_N50_biotin_gene))
print(up_shared)

### Comparing expression levels across Input and Biotin samples
# N50
Downregulated_N50_IP <- N50_input_df %>% filter(Expression == "Downregulated") %>% select("logFC", "gene_id")
colnames(Downregulated_N50_IP)[1] <- "IP_N50"
Downregulated_N50_Biotin <- N50_biotin_df %>% filter(Expression == "Downregulated") %>% select("logFC", "gene_id")
colnames(Downregulated_N50_Biotin)[1] <- "Biotin_N50"
join_N50 <- inner_join(Downregulated_N50_IP, Downregulated_N50_Biotin, by="gene_id")

# Creating the plot with positive values
plot(abs(join_N50$Biotin_N50), abs(join_N50$IP_N50), 
     xlab = "abs(log2(Tbp1∆50 Biotin / WT)", ylab = "abs(log2(Tbp1∆50 Input / WT)",
     main = "Comparison of logFC between Input and Biotin Samples",
     pch = 16, col = "black")  # Customize plot symbols and color if needed

# Regression line
abline(lm(abs(join_N50$IP_N50) ~ abs(join_N50$Biotin_N50)), col = "red", lwd = 3)

# Compute Pearson correlation coefficient
correlation <- cor(abs(join_N50$IP_N50), abs(join_N50$Biotin_N50))

# Add Pearson correlation as text
text(3, 3, paste("Correlation:", round(correlation, 2)), col = "red")


# N25
Downregulated_N25_IP <- N25_input_df %>% filter(Expression == "Downregulated") %>% select("logFC", "gene_id")
colnames(Downregulated_N25_IP)[1] <- "IP_N25"
Downregulated_N25_Biotin <- N25_biotin_df %>% filter(Expression == "Downregulated") %>% select("logFC", "gene_id")
colnames(Downregulated_N25_Biotin)[1] <- "Biotin_N25"
join_N25 <- inner_join(Downregulated_N25_IP, Downregulated_N25_Biotin, by="gene_id")

# Creating the plot with positive values
plot(abs(join_N25$Biotin_N25), abs(join_N25$IP_N25), 
     xlab = "abs(log2(Tbp1∆50 Biotin / WT))", ylab = "abs(log2(Tbp1∆50 Input / WT))",
     main = "Comparison of logFC between Input and Biotin Samples",
     pch = 16, col = "black")  # Customize plot symbols and color if needed

# Regression line
abline(lm(abs(join_N25$IP_N25) ~ abs(join_N25$Biotin_N25)), col = "red", lwd = 3)

# Compute Pearson correlation coefficient
correlation <- cor(abs(join_N25$IP_N25), abs(join_N25$Biotin_N25))

# Add Pearson correlation as text
text(1.5, 1.5, paste("Correlation:", round(correlation, 2)), col = "red")


