# Set working directory
setwd("Data/")

# Load packages and Seurat objects
library(Seurat)
library(dplyr)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)

xing <- readRDS("Xing2020/xing_final.rds")
pana <- readRDS("Panariello2023/pana_merged_final.rds")


# Identify cells that express high CDKN1A or CDKN2A
exp_cdkn1a.xing <- WhichCells(xing, expression = CDKN1A > 1)
exp_cdkn2a.xing <- WhichCells(xing, expression = CDKN2A > 1)
exp_sen.xing <- union(exp_cdkn1a.xing, exp_cdkn2a.xing)

xing$ExpSenGenes <- NA
xing$ExpSenGenes <- ifelse(colnames(xing) %in% exp_sen.xing, "Present", "Absent")

exp_cdkn1a.pana <- WhichCells(pana, expression = CDKN1A > 1)
exp_cdkn2a.pana <- WhichCells(pana, expression = CDKN2A > 1)
exp_sen.pana <- union(exp_cdkn1a.pana, exp_cdkn2a.pana)

pana$ExpSenGenes <- NA
pana$ExpSenGenes <- ifelse(colnames(pana) %in% exp_sen.pana, "Present", "Absent")


# Collapse single-cell data in pseudo-bulk
pb.xing <- AggregateExpression(
  object = xing, 
  group.by = "ExpSenGenes", 
  assays = "RNA", 
  slot = "counts"
)
pb.xing <- pb.xing[[1]]

pb.pana <- AggregateExpression(
  object = pana, 
  group.by = "ExpSenGenes", 
  assays = "RNA", 
  slot = "counts"
)
pb.pana <- pb.pana[[1]]


## Normalize pseudobulk counts with and perform DE analysis w/ DESeq2
# Create matrix
genes <- union(rownames(pb.xing), rownames(pb.pana))
counts <- matrix(
  data = 0,
  nrow = length(genes), ncol = 4, 
  dimnames = list(genes, c("Xing_Pro", "Xing_Sen", "Panariello_Pro", "Panariello_Sen"))
)

for (gene in rownames(counts)) {
  if (gene %in% rownames(pb.xing)) {
    counts[gene, "Xing_Pro"] <- pb.xing[gene, "Absent"]
    counts[gene, "Xing_Sen"] <- pb.xing[gene, "Present"]
  }
  
  if (gene %in% rownames(pb.pana)) {
    counts[gene, "Panariello_Pro"] <- pb.pana[gene, "Absent"]
    counts[gene, "Panariello_Sen"] <- pb.pana[gene, "Present"]
  }
}

metadata <- data.frame(
  sample = c("Xing", "Xing", "Panariello", "Panariello"),
  condition = c("Proliferating", "Senescent", "Proliferating", "Senescent"), 
  row.names = colnames(counts)
)

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ sample + condition
)
dds <- DESeq(dds)
norm.counts <- counts(dds, normalized = TRUE)
res <- lfcShrink(
  dds = dds,
  coef = "condition_Senescent_vs_Proliferating", 
  type = "apeglm"
)
res <- res[complete.cases(res),]
res <- res[order(res$log2FoldChange, decreasing = TRUE),]

# Save table of all significantly 
# differentially expressed genes in RIS, 
# ranked by log-fold change
write.table(
  as.data.frame(res), 
  file = "../Tables/SenotypeComparisons/sen_all_normalized_repro.tsv", 
  sep = "\t", qmethod = "double"
)


## Save list of genes used in DESeq for background gene set during GO/KEGG analysis
# Convert gene names to Entrez IDs
bg.genes <- bitr(
  rownames(res), 
  fromType = c("SYMBOL"),
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db, 
  drop = FALSE
)

# Manually input Entrez IDs of RIS marker genes that failed to map due to alternate names
bg.genes[which(bg.genes$SYMBOL == "LINC01137"), "ENTREZID"] <- 728431
bg.genes[which(bg.genes$SYMBOL == "KIAA1324"), "ENTREZID"] <- 57535
bg.genes[which(bg.genes$SYMBOL == "AHSA2"), "ENTREZID"] <- 130872
bg.genes[which(bg.genes$SYMBOL == "FAM26E"), "ENTREZID"] <- 254228
bg.genes[which(bg.genes$SYMBOL == "FAM188B"), "ENTREZID"] <- 84182
bg.genes[which(bg.genes$SYMBOL == "C11orf70"), "ENTREZID"] <- 85016
bg.genes[which(bg.genes$SYMBOL == "LYRM5"), "ENTREZID"] <- 144363
bg.genes[which(bg.genes$SYMBOL == "METTL7A"), "ENTREZID"] <- 25840
bg.genes[which(bg.genes$SYMBOL == "LHFP"), "ENTREZID"] <- 10186
bg.genes[which(bg.genes$SYMBOL == "C16orf45"), "ENTREZID"] <- 89927
bg.genes[which(bg.genes$SYMBOL == "C17orf82"), "ENTREZID"] <- 388407
bg.genes[which(bg.genes$SYMBOL == "MARCH2"), "ENTREZID"] <- 51257
bg.genes[which(bg.genes$SYMBOL == "WDR63"), "ENTREZID"] <- 126820
bg.genes[which(bg.genes$SYMBOL == "HIST2H2BE"), "ENTREZID"] <- 8349
bg.genes[which(bg.genes$SYMBOL == "FAM198B"), "ENTREZID"] <- 51313
bg.genes[which(bg.genes$SYMBOL == "C5orf56"), "ENTREZID"] <- 441108
bg.genes[which(bg.genes$SYMBOL == "FAM65B"), "ENTREZID"] <- 9750
bg.genes[which(bg.genes$SYMBOL == "HIST1H2AC"), "ENTREZID"] <- 8334
bg.genes[which(bg.genes$SYMBOL == "HIST1H2BD"), "ENTREZID"] <- 3017
bg.genes[which(bg.genes$SYMBOL == "C8orf4"), "ENTREZID"] <- 56892
bg.genes[which(bg.genes$SYMBOL == "FAM214B"), "ENTREZID"] <- 80256
bg.genes[which(bg.genes$SYMBOL == "OBFC1"), "ENTREZID"] <- 79991
bg.genes[which(bg.genes$SYMBOL == "RARRES3"), "ENTREZID"] <- 5920
bg.genes[which(bg.genes$SYMBOL == "LINC00346"), "ENTREZID"] <- 283487
bg.genes[which(bg.genes$SYMBOL == "FAM214A"), "ENTREZID"] <- 56204
bg.genes[which(bg.genes$SYMBOL == "FAM65C"), "ENTREZID"] <- 140876
bg.genes[which(bg.genes$SYMBOL == "RP3-510D11.2"), "ENTREZID"] <- 102724571
bg.genes[which(bg.genes$SYMBOL == "RP11-14N7.2"), "ENTREZID"] <- 100996732
bg.genes[which(bg.genes$SYMBOL == "RP11-400N13.3"), "ENTREZID"] <- 105372950
bg.genes[which(bg.genes$SYMBOL == "AC007246.3"), "ENTREZID"] <- 728730
bg.genes[which(bg.genes$SYMBOL == "AC159540.1"), "ENTREZID"] <- 100506076
bg.genes[which(bg.genes$SYMBOL == "RP4-613B23.1"), "ENTREZID"] <- 101928323
bg.genes[which(bg.genes$SYMBOL == "RP11-539L10.3"), "ENTREZID"] <- 105374366
bg.genes[which(bg.genes$SYMBOL == "RP11-597D13.9"), "ENTREZID"] <- 285505
bg.genes[which(bg.genes$SYMBOL == "RP11-472N13.3"), "ENTREZID"] <- 101929951
bg.genes[which(bg.genes$SYMBOL == "RP11-783K16.5"), "ENTREZID"] <- 105369340
bg.genes[which(bg.genes$SYMBOL == "RP11-658F2.8"), "ENTREZID"] <- 5936
bg.genes[which(bg.genes$SYMBOL == "RP11-362K2.2"), "ENTREZID"] <- 100506869
bg.genes[which(bg.genes$SYMBOL == "RP11-66N24.3"), "ENTREZID"] <- 102724814
bg.genes[which(bg.genes$SYMBOL == "AC006262.4"), "ENTREZID"] <- 105372424
bg.genes[which(bg.genes$SYMBOL == "AC006262.5"), "ENTREZID"] <- 645553

# Remove genes that did not have Entrez IDs
bg.genes <- bg.genes[which(!is.na(bg.genes$ENTREZID)),]

print(paste0(
  nrow(bg.genes), "/", 
  nrow(res), 
  " genes used in the DESeq testing had Entrez IDs."
))

write.table(
  bg.genes, file = "../Tables/SenotypeComparisons/repro_geneIDs_background.tsv", 
  quote = 1, sep = "\t", 
  row.names = FALSE, qmethod = "double"
)


## Identify genes upregulated and downregulated in senescence
res <- res[which(res$padj < 0.05),]
up.genes <- rownames(res[which(res$log2FoldChange > 0.5),])
up.counts <- norm.counts[which(rownames(norm.counts) %in% up.genes),]
down.genes <- rownames(res[which(res$log2FoldChange < -0.5),])
down.counts <- norm.counts[which(rownames(norm.counts) %in% down.genes),]

# Save table of significantly upregulated genes during senescence
write.table(
  up.counts, file = "../Tables/SenotypeComparisons/sen_upgenes_normalized_repro.tsv", 
  sep = "\t", qmethod = "double"
)
write.table(
  down.counts, file = "../Tables/SenotypeComparisons/sen_downgenes_normalized_repro.tsv", 
  sep = "\t", qmethod = "double"
)


