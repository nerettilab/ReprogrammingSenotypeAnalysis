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
res <- results(dds)


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
bg.genes[which(bg.genes$SYMBOL == "LINC00152"), "ENTREZID"] <- 112597
bg.genes[which(bg.genes$DYMBOL == "GPR1"), "ENTREZID"] <- 2825
bg.genes[which(bg.genes$SYMBOL == "FLJ43879"), "ENTREZID"] <- 401039
bg.genes[which(bg.genes$SYMBOL == "KIAA1109"), "ENTREZID"] <- 84162
bg.genes[which(bg.genes$SYMBOL == "SEPT8"), "ENTREZID"] <- 23176
bg.genes[which(bg.genes$SYMBOL == "FAM26E"), "ENTREZID"] <- 254228
bg.genes[which(bg.genes$SYMBOL == "FAM188B"), "ENTREZID"] <- 84182
bg.genes[which(bg.genes$SYMBOL == "ABHD11-AS1"), "ENTREZID"] <- 171022
bg.genes[which(bg.genes$SYMBOL == "FAM71F2"), "ENTREZID"] <- 346653
bg.genes[which(bg.genes$SYMBOL == "CMB9-22P13.1"), "ENTREZID"] <- 101927789
bg.genes[which(bg.genes$SYMBOL == "C11orf70"), "ENTREZID"] <- 85016
bg.genes[which(bg.genes$SYMBOL == "LYRM5"), "ENTREZID"] <- 144363
bg.genes[which(bg.genes$SYMBOL == "METTL7A"), "ENTREZID"] <- 25840
bg.genes[which(bg.genes$SYMBOL == "LHFP"), "ENTREZID"] <- 10186
bg.genes[which(bg.genes$SYMBOL == "C15orf52"), "ENTREZID"] <- 388115
bg.genes[which(bg.genes$SYMBOL == "C16orf45"), "ENTREZID"] <- 89927
bg.genes[which(bg.genes$SYMBOL == "FAM134C"), "ENTREZID"] <- 162427
bg.genes[which(bg.genes$SYMBOL == "C17orf82"), "ENTREZID"] <- 388407
bg.genes[which(bg.genes$SYMBOL == "BAIAP2-AS1"), "ENTREZID"] <- 440465
bg.genes[which(bg.genes$SYMBOL == "ZADH2"), "ENTREZID"] <- 284273
bg.genes[which(bg.genes$SYMBOL == "MMP24-AS1"), "ENTREZID"] <- 101410538
bg.genes[which(bg.genes$SYMBOL == "ATP5E"), "ENTREZID"] <- 514
bg.genes[which(bg.genes$SYMBOL == "MARCH2"), "ENTREZID"] <- 51257
bg.genes[which(bg.genes$SYMBOL == "WDR63"), "ENTREZID"] <- 126820
bg.genes[which(bg.genes$SYMBOL == "HIST2H2BE"), "ENTREZID"] <- 8349
bg.genes[which(bg.genes$SYMBOL == "SDPR"), "ENTREZID"] <- 8436
bg.genes[which(bg.genes$SYMBOL == "FAM198B"), "ENTREZID"] <- 51313
bg.genes[which(bg.genes$SYMBOL == "C5orf56"), "ENTREZID"] <- 441108
bg.genes[which(bg.genes$SYMBOL == "FAM65B"), "ENTREZID"] <- 9750
bg.genes[which(bg.genes$SYMBOL == "HIST1H1C"), "ENTREZID"] <- 3006
bg.genes[which(bg.genes$SYMBOL == "HIST1H2AC"), "ENTREZID"] <- 8334
bg.genes[which(bg.genes$SYMBOL == "HIST1H2BD"), "ENTREZID"] <- 3017
bg.genes[which(bg.genes$SYMBOL == "C8orf4"), "ENTREZID"] <- 56892
bg.genes[which(bg.genes$SYMBOL == "FAM214B"), "ENTREZID"] <- 80256
bg.genes[which(bg.genes$SYMBOL == "C10orf54"), "ENTREZID"] <- 64115
bg.genes[which(bg.genes$SYMBOL == "OBFC1"), "ENTREZID"] <- 79991
bg.genes[which(bg.genes$SYMBOL == "RARRES3"), "ENTREZID"] <- 5920
bg.genes[which(bg.genes$SYMBOL == "LINC00346"), "ENTREZID"] <- 283487
bg.genes[which(bg.genes$SYMBOL == "FAM214A"), "ENTREZID"] <- 56204
bg.genes[which(bg.genes$SYMBOL == "FAM65C"), "ENTREZID"] <- 140876

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


## Identify markers of senescence
res <- res[which(res$padj < 0.05),]
top.genes <- rownames(res[which(res$log2FoldChange > 0.5),])
top.counts <- norm.counts[which(rownames(norm.counts) %in% top.genes),]

# Save table of significantly upregulated genes during senescence
write.table(
  top.counts, file = "../Tables/SenotypeComparisons/sen_genes_normalized_repro.tsv", 
  sep = "\t", qmethod = "double"
)


