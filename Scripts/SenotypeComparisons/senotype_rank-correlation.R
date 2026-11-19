# Set working directory
setwd("Tables/SenotypeComparisons/")

# Load packages and counts tables
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)

bj.repro <- read.delim("sen_all_normalized_repro.tsv")

bj.rs <- read.csv(file = "../../Data/SenotypeComparisons/BJ_RS_marthandan2016_counts.csv")

bj.sis.young.9_11 <- read.table(file = "../../Data/SenotypeComparisons/BJ_SIS_alspach2014_young_9_11.txt", sep = "\t", col.names = c("external_gene_id", "count"))
bj.sis.young.9_24 <- read.table(file = "../../Data/SenotypeComparisons/BJ_SIS_alspach2014_young_9_24.txt", sep = "\t", col.names = c("external_gene_id", "count"))
bj.sis.sen.9_11 <- read.table(file = "../../Data/SenotypeComparisons/BJ_SIS_alspach2014_sen_9_11.txt", sep = "\t", col.names = c("external_gene_id", "count"))
bj.sis.sen.9_24 <- read.table(file = "../../Data/SenotypeComparisons/BJ_SIS_alspach2014_sen_9_24.txt", sep = "\t", col.names = c("external_gene_id", "count"))

bj.ois <- read.delim("../../Data/SenotypeComparisons/BJ_OIS_loayza-puch2013_counts.tsv")
colnames(bj.ois) <- c(
  "GeneID", "Control_1", "Control_2", "Quiesc.SW_1", 
  "Quiesc.SW_2", "RASG12V.Day05_1", "RASG12V.Day05_2", 
  "Control_3", "RASG12V.Day14", "Transformed"
)


## Normalize counts and perform DE analysis w/ DESeq2
# BJ - RS
# Manually rename non-unique IDs for genes that will end up in results data frame
bj.rs[bj.rs$ensembl_gene_id == "ENSG00000254154", "external_gene_id"] <- "CRYZL2P-SEC16B"
bj.rs[bj.rs$ensembl_gene_id == "ENSG00000173867", "external_gene_id"] <- "ENSG00000173867"

counts.bj.rs <- bj.rs[,c("BJ_Y1", "BJ_Y2", "BJ_Y3", "BJ_OLD_1", "BJ_OLD_2", "BJ_OLD_3")]
rownames(counts.bj.rs) <- bj.rs$ensembl_gene_id
metadata.bj.rs <- data.frame(
  condition = c("Proliferating", "Proliferating", "Proliferating", "Senescent", "Senescent", "Senescent"), 
  row.names = colnames(counts.bj.rs)
)

dds.bj.rs <- DESeqDataSetFromMatrix(
  countData = counts.bj.rs, 
  colData = metadata.bj.rs, 
  design = ~ condition
)
dds.bj.rs <- DESeq(dds.bj.rs)
norm.counts.bj.rs <- counts(dds.bj.rs, normalized = TRUE)
norm.counts.bj.rs <- cbind(
  norm.counts.bj.rs, 
  external_gene_id = bj.rs$external_gene_id
)

res.bj.rs <- lfcShrink(
  dds = dds.bj.rs,
  coef = "condition_Senescent_vs_Proliferating"
)
res.bj.rs <- as.data.frame(res.bj.rs[order(res.bj.rs$log2FoldChange, decreasing = TRUE),])
res.bj.rs$external_gene_id <- NA
for (row in rownames(res.bj.rs)) {
  res.bj.rs[row, "external_gene_id"] <- bj.rs[bj.rs$ensembl_gene_id == row, "external_gene_id"]
}
res.bj.rs <- res.bj.rs[res.bj.rs$external_gene_id != "",]
res.bj.rs <- res.bj.rs[complete.cases(res.bj.rs),]
norm.counts.bj.rs <- as.data.frame(norm.counts.bj.rs[which(rownames(norm.counts.bj.rs) %in% rownames(res.bj.rs)),])

# BJ - SIS
counts.bj.sis <- matrix(
  data = c(bj.sis.young.9_11$count, bj.sis.young.9_24$count, bj.sis.sen.9_11$count, bj.sis.sen.9_24$count), 
  nrow = nrow(bj.sis.young.9_11), ncol = 4, 
  dimnames = list(bj.sis.young.9_11$external_gene_id, c("Young.9_11", "Young.9_24", "Sen.9_11", "Sen.9_24"))
)
metadata.bj.sis <- data.frame(
  condition = c("Proliferating", "Proliferating", "Senescent", "Senescent"), 
  row.names = colnames(counts.bj.sis)
)

dds.bj.sis <- DESeqDataSetFromMatrix(
  countData = counts.bj.sis,
  colData = metadata.bj.sis,
  design = ~ condition
)
dds.bj.sis <- DESeq(dds.bj.sis)
norm.counts.bj.sis <- counts(dds.bj.sis, normalized = TRUE)

res.bj.sis <- lfcShrink(
  dds = dds.bj.sis,
  coef = "condition_Senescent_vs_Proliferating", 
  type = "apeglm"
)
res.bj.sis <- as.data.frame(res.bj.sis[order(res.bj.sis$log2FoldChange, decreasing = TRUE),])
res.bj.sis <- res.bj.sis[complete.cases(res.bj.sis),]
norm.counts.bj.sis <- as.data.frame(norm.counts.bj.sis[which(rownames(norm.counts.bj.sis) %in% rownames(res.bj.sis)),])

# BJ - OIS
counts.bj.ois <- matrix(
  data = c(bj.ois$Control_1, bj.ois$Control_2, bj.ois$Control_3, bj.ois$RASG12V.Day14), 
  nrow = nrow(bj.ois), ncol = 4, 
  dimnames = list(bj.ois$GeneID, c("Control_1", "Control_2", "Control_3", "Senescent"))
)
metadata.bj.ois <- data.frame(
  condition = c("Proliferating", "Proliferating", "Proliferating", "Senescent"), 
  row.names = colnames(counts.bj.ois)
)

dds.bj.ois <- DESeqDataSetFromMatrix(
  countData = counts.bj.ois,
  colData = metadata.bj.ois,
  design = ~ condition
)
dds.bj.ois <- DESeq(dds.bj.ois)
norm.counts.bj.ois <- counts(dds.bj.ois, normalized = TRUE)
gene.symbols.bj.ois <- mapIds(
  x = org.Hs.eg.db, keys = as.character(bj.ois$GeneID), 
  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
)
norm.counts.bj.ois <- cbind(
  norm.counts.bj.ois, 
  GeneID = gene.symbols.bj.ois
)

res.bj.ois <- lfcShrink(
  dds = dds.bj.ois,
  coef = "condition_Senescent_vs_Proliferating", 
  type = "apeglm"
)
res.bj.ois <- as.data.frame(res.bj.ois[order(res.bj.ois$log2FoldChange, decreasing = TRUE),])
res.bj.ois$GeneID <- NA
for (row in rownames(res.bj.ois)) {
  res.bj.ois[row, "GeneID"] <- norm.counts.bj.ois[rownames(norm.counts.bj.ois) == row, "GeneID"]
}
res.bj.ois <- res.bj.ois[complete.cases(res.bj.ois),]
norm.counts.bj.ois <- as.data.frame(norm.counts.bj.ois[which(rownames(norm.counts.bj.ois) %in% rownames(res.bj.ois)),])


## Rank-correlate log-fold change vectors 
## across senotypes
all.comm <- Reduce(
  intersect, list(
    rownames(bj.repro),
    unique(res.bj.rs$external_gene_id), 
    rownames(res.bj.sis), 
    unique(res.bj.ois$GeneID)
  )
)
lfc.ranks <- data.frame(matrix(
  ncol = 4, nrow = length(all.comm), 
  dimnames = list(
    all.comm, 
    c("Reprogramming", "Replicative", "Stress", "Oncogene")
  )
))

for (row in rownames(lfc.ranks)) {
  lfc.ranks[row, "Reprogramming"] <- bj.repro[row, "log2FoldChange"]
  lfc.ranks[row, "Replicative"] <- res.bj.rs[res.bj.rs$external_gene_id == row, "log2FoldChange"]
  lfc.ranks[row, "Stress"] <- res.bj.sis[row, "log2FoldChange"]
  lfc.ranks[row, "Oncogene"] <- res.bj.ois[res.bj.ois$GeneID == row, "log2FoldChange"]
}

lfc.corr <- data.frame(matrix(
  ncol = 4, nrow = 4, 
  dimnames = list(colnames(lfc.ranks), colnames(lfc.ranks))
))
padj.corr <- lfc.corr

comps <- combn(colnames(lfc.ranks), m = 2)
for (i in 1:ncol(comps)) {
  sen1 <- comps[1,i]
  sen2 <- comps[2,i]
  
  corr <- cor.test(
    x = lfc.ranks[,sen1], 
    y = lfc.ranks[,sen2],
    method = "spearman"
  )
  
  lfc.corr[sen1, sen2] <- corr$estimate
  padj.corr[sen1, sen2] <- corr$p.value
} # Create tables of Spearman coefficients and p-values

padj <- p.adjust(padj.corr[!is.na(padj.corr)], method = "bonferroni")
for (i in 1:ncol(comps)) {
  sen1 <- comps[1,i]
  sen2 <- comps[2,i]
  
  padj.corr[sen1, sen2] <- padj[i]
} # Bonferroni correction of p-values

lfc.corr

# Robustness check: rank-correlation after 
# filtering out genes not appreciably 
# differentially expressed (|log2FC| <= 0.5)
# in at least half of the senotype
keep <- rowMeans(abs(lfc.ranks) > 0.5) >= 0.5
lfc.corr.exp <- data.frame(matrix(
  ncol = 4, nrow = 4, 
  dimnames = list(colnames(lfc.ranks), colnames(lfc.ranks))
))
padj.corr.exp <- lfc.corr.exp

for (i in 1:ncol(comps)) {
  sen1 <- comps[1,i]
  sen2 <- comps[2,i]
  
  corr <- cor.test(
    x = lfc.ranks[keep, sen1], 
    y = lfc.ranks[keep, sen2],
    method = "spearman"
  )
  
  lfc.corr.exp[sen1, sen2] <- corr$estimate
  padj.corr.exp[sen1, sen2] <- corr$p.value
} # Create tables of Spearman coefficients and p-values

padj <- p.adjust(padj.corr.exp[!is.na(padj.corr.exp)], method = "bonferroni")
for (i in 1:ncol(comps)) {
  sen1 <- comps[1,i]
  sen2 <- comps[2,i]
  
  padj.corr.exp[sen1, sen2] <- padj[i]
} # Bonferroni correction of p-values

lfc.corr.exp

# Robustness check: rank-correlation after 
# filtering out genes whose log2-fold changes 
# don't have p-values below 0.05 in at least 
# half of the senotypes
sig.genes <- vector(mode = "logical", length = length(all.comm))
names(sig.genes) <- all.comm
for (gene in names(sig.genes)) {
  pvals <- c(
    bj.repro[gene, "padj"], 
    res.bj.rs[res.bj.rs$external_gene_id == gene, "padj"], 
    res.bj.sis[gene, "padj"], 
    res.bj.ois[res.bj.ois$GeneID == row, "padj"]
  )
  
  sig.genes[gene] <- mean(pvals < 0.05) == 0.5
}

lfc.corr.sig <- data.frame(matrix(
  ncol = 4, nrow = 4, 
  dimnames = list(colnames(lfc.ranks), colnames(lfc.ranks))
))
padj.corr.sig <- lfc.corr.sig

for (i in 1:ncol(comps)) {
  sen1 <- comps[1,i]
  sen2 <- comps[2,i]
  
  corr <- cor.test(
    x = lfc.ranks[sig.genes, sen1], 
    y = lfc.ranks[sig.genes, sen2],
    method = "spearman"
  )
  
  lfc.corr.sig[sen1, sen2] <- corr$estimate
  padj.corr.sig[sen1, sen2] <- corr$p.value
} # Create tables of Spearman coefficients and p-values

padj <- p.adjust(padj.corr.sig[!is.na(padj.corr.sig)], method = "bonferroni")
for (i in 1:ncol(comps)) {
  sen1 <- comps[1,i]
  sen2 <- comps[2,i]
  
  padj.corr.sig[sen1, sen2] <- padj[i]
} # Bonferroni correction of p-values

lfc.corr.sig


