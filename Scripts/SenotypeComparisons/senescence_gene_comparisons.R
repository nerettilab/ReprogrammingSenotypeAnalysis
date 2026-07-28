# Set working directory
setwd("Tables/SenotypeComparisons/")

# Load packages and counts tables
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(ggupset)

bj.repro <- read.delim("sen_genes_normalized_repro.tsv")

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

res.bj.rs <- results(dds.bj.rs)
res.bj.rs <- res.bj.rs[which(res.bj.rs$padj < 0.05),]
top.genes.rs <- rownames(res.bj.rs[which(res.bj.rs$log2FoldChange > 0.5),])
top.bj.rs <- norm.counts.bj.rs[which(rownames(norm.counts.bj.rs) %in% top.genes.rs),]
top.bj.rs <- top.bj.rs[which(top.bj.rs[,"external_gene_id"] != ""),]

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

res.bj.sis <- results(dds.bj.sis)
res.bj.sis <- res.bj.sis[which(res.bj.sis$padj < 0.05),]
top.genes.sis <- rownames(res.bj.sis[which(res.bj.sis$log2FoldChange > 0.5),])
top.bj.sis <- norm.counts.bj.sis[which(rownames(norm.counts.bj.sis) %in% top.genes.sis),]

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

res.bj.ois <- results(dds.bj.ois)
res.bj.ois <- res.bj.ois[which(res.bj.ois$padj < 0.05),]
top.genes.ois <- rownames(res.bj.ois[which(res.bj.ois$log2FoldChange > 0.5),])
top.bj.ois <- norm.counts.bj.ois[which(rownames(norm.counts.bj.ois) %in% top.genes.ois),]
top.bj.ois <- top.bj.ois[which(top.bj.ois[,"GeneID"] != ""),]


## Compare markers of reprogramming-induced senescent cells with those of other senotypes
# BJ cells (Human)
repro.rs.bj <- which(rownames(bj.repro) %in% top.bj.rs[,"external_gene_id"])
repro.sis.bj <- which(rownames(bj.repro) %in% rownames(top.bj.sis))
repro.ois.bj <- which(rownames(bj.repro) %in% top.bj.ois[,"GeneID"])

print(paste0(
  length(union(repro.rs.bj, union(repro.sis.bj, repro.ois.bj))), "/", 
  length(rownames(bj.repro)), 
  " markers for reprogramming-induced senescence are also among the markers for another type of senescence induction."
))
print(paste0(
  length(rownames(bj.repro)) - length(union(repro.rs.bj, union(repro.sis.bj, repro.ois.bj))), "/", 
  length(rownames(bj.repro)), 
  " markers for reprogramming-induced senescence appear to be unique to this form of senescence induction."
))

print(paste(
  length(repro.rs.bj), 
  "markers for reprogramming-induced senescence are found among the markers for RS."
))
print(paste(
  length(repro.sis.bj), 
  "markers for reprogramming-induced senescence are found among the markers for SIS."
))
print(paste(
  length(repro.ois.bj), 
  "markers for reprogramming-induced senescence are found among the markers for OIS."
))

print(paste(
  length(intersect(repro.rs.bj, repro.sis.bj)), 
  "markers for reprogramming-induced senescence are found among markers common to RS and SIS."
))
print(paste(
  length(intersect(repro.rs.bj, repro.ois.bj)), 
  "markers for reprogramming-induced senescence are found among markers common to RS and OIS."
))
print(paste(
  length(intersect(repro.sis.bj, repro.ois.bj)), 
  "markers for reprogramming-induced senescence are found among markers common to SIS and OIS."
))
print(paste(
  length(Reduce(intersect, list(repro.rs.bj, repro.sis.bj, repro.ois.bj))), 
  "markers for reprogramming-induced senescence are also common to RS, SIS, and OIS."
))


## UpSet plot for BJ cells
# Create tables of reprogramming-induced SenSig markers and their presence as markers of other senescence types
all.sen.genes <- Reduce(union, list(
  rownames(bj.repro), 
  unique(top.bj.rs[,"external_gene_id"]), 
  rownames(top.bj.sis), 
  unique(top.bj.ois[,"GeneID"])
))
bj.senmarkers <- data.frame(matrix(
  data = FALSE, 
  nrow = length(all.sen.genes), 
  ncol = 5
))
colnames(bj.senmarkers) <- c("Gene", "Repro", "RS", "SIS", "OIS")
bj.senmarkers$Gene <- all.sen.genes

for (gene in 1:nrow(bj.senmarkers)) {
  gene.name <- bj.senmarkers$Gene[gene]
  
  bj.senmarkers[gene, "Repro"] <- ifelse(
    gene.name %in% rownames(bj.repro), 
    TRUE, FALSE
  )
  bj.senmarkers[gene, "RS"] <- ifelse(
    gene.name %in% unique(top.bj.rs[,"external_gene_id"]), 
    TRUE, FALSE
  )
  bj.senmarkers[gene, "SIS"] <- ifelse(
    gene.name %in% rownames(top.bj.sis), 
    TRUE, FALSE
  )
  bj.senmarkers[gene, "OIS"] <- ifelse(
    gene.name %in% unique(top.bj.ois[,"GeneID"]), 
    TRUE, FALSE
  )
}

# Format table for use in ggupset
tidy.bj.senmarkers <- as_tibble(bj.senmarkers)
tidy.bj.senmarkers <- tidy.bj.senmarkers |> mutate(
  combination = pmap(
    list(Repro, RS, SIS, OIS), 
    \(lgl1, lgl2, lgl3, lgl4) {
      c("Reprogramming", "Replicative", "Stress", "Oncogene")[c(lgl1, lgl2, lgl3, lgl4)]
    }
  )
)

# Create and save UpSet plots
bj.us <- tidy.bj.senmarkers |> 
  ggplot(aes(x = combination)) + geom_bar() + geom_text(
    stat = "count", aes(label = after_stat(count)), 
    vjust = -1, size = 6
  ) + scale_x_upset() + labs(
    x = element_blank(), 
    y = "# of Genes"
  ) + theme(
    legend.position = "none", 
    text = element_text(size = 25)
  ) + theme_combmatrix(
    combmatrix.label.extra_spacing = 1, 
    combmatrix.label.text = element_text(size = 15), 
    combmatrix.panel.point.size = 5, 
    combmatrix.panel.line.size = 1.7
  )

png(
  filename = "../../Images/SenotypeComparisons/upsetplot_senmarkers.png", 
  width = 1680, height = 1800, res = 200
)
bj.us
dev.off()

# Save table of marker genes unique to reprogramming-induced senescence
ris.specific <- vector(mode = "numeric", length = 439)
for (gene in 1:nrow(tidy.bj.senmarkers)) {
  if (identical(tidy.bj.senmarkers$combination[[gene]], c("Reprogramming"))) {
    ris.specific[match(0, ris.specific)] <- gene
  }
}

specific.genes <- as.data.frame(tidy.bj.senmarkers[ris.specific, -ncol(tidy.bj.senmarkers)])
write.table(
  specific.genes, 
  file = "senotype_reprogramming-specific_genes.tsv", 
  row.names = FALSE, sep = "\t", qmethod = "double"
)

# Save table of marker genes common to all 4 forms of senescence induction
all4 <- vector(mode = "numeric", length = 32)
for (gene in 1:nrow(tidy.bj.senmarkers)) {
  if (identical(tidy.bj.senmarkers$combination[[gene]], c("Reprogramming", "Replicative", "Stress", "Oncogene"))) {
    all4[match(0, all4)] <- gene
  }
}

common.genes <- as.data.frame(tidy.bj.senmarkers[all4,-ncol(tidy.bj.senmarkers)])
write.table(
  common.genes, 
  file = "senotype_common_genes.tsv", 
  row.names = FALSE, sep = "\t", qmethod = "double"
)


## GO/KEGG Analysis of RIS-specific marker genes
## (All RIS markers)
# Create new gene set
genes.df <- bitr(
  rownames(bj.repro),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = FALSE
)

# Manually input Entrez IDs of genes that failed to map due to alternate names
genes.df[which(genes.df$SYMBOL == "LINC01137"), "ENTREZID"] <- 728431
genes.df[which(genes.df$SYMBOL == "KIAA1324"), "ENTREZID"] <- 57535
genes.df[which(genes.df$SYMBOL == "AHSA2"), "ENTREZID"] <- 130872
genes.df[which(genes.df$SYMBOL == "LINC00152"), "ENTREZID"] <- 112597
genes.df[which(genes.df$DYMBOL == "GPR1"), "ENTREZID"] <- 2825
genes.df[which(genes.df$SYMBOL == "FLJ43879"), "ENTREZID"] <- 401039
genes.df[which(genes.df$SYMBOL == "KIAA1109"), "ENTREZID"] <- 84162
genes.df[which(genes.df$SYMBOL == "SEPT8"), "ENTREZID"] <- 23176
genes.df[which(genes.df$SYMBOL == "FAM26E"), "ENTREZID"] <- 254228
genes.df[which(genes.df$SYMBOL == "FAM188B"), "ENTREZID"] <- 84182
genes.df[which(genes.df$SYMBOL == "ABHD11-AS1"), "ENTREZID"] <- 171022
genes.df[which(genes.df$SYMBOL == "FAM71F2"), "ENTREZID"] <- 346653
genes.df[which(genes.df$SYMBOL == "CMB9-22P13.1"), "ENTREZID"] <- 101927789
genes.df[which(genes.df$SYMBOL == "C11orf70"), "ENTREZID"] <- 85016
genes.df[which(genes.df$SYMBOL == "LYRM5"), "ENTREZID"] <- 144363
genes.df[which(genes.df$SYMBOL == "METTL7A"), "ENTREZID"] <- 25840
genes.df[which(genes.df$SYMBOL == "LHFP"), "ENTREZID"] <- 10186
genes.df[which(genes.df$SYMBOL == "C15orf52"), "ENTREZID"] <- 388115
genes.df[which(genes.df$SYMBOL == "C16orf45"), "ENTREZID"] <- 89927
genes.df[which(genes.df$SYMBOL == "FAM134C"), "ENTREZID"] <- 162427
genes.df[which(genes.df$SYMBOL == "C17orf82"), "ENTREZID"] <- 388407
genes.df[which(genes.df$SYMBOL == "BAIAP2-AS1"), "ENTREZID"] <- 440465
genes.df[which(genes.df$SYMBOL == "ZADH2"), "ENTREZID"] <- 284273
genes.df[which(genes.df$SYMBOL == "MMP24-AS1"), "ENTREZID"] <- 101410538
genes.df[which(genes.df$SYMBOL == "ATP5E"), "ENTREZID"] <- 514
genes.df[which(genes.df$SYMBOL == "MARCH2"), "ENTREZID"] <- 51257
genes.df[which(genes.df$SYMBOL == "WDR63"), "ENTREZID"] <- 126820
genes.df[which(genes.df$SYMBOL == "HIST2H2BE"), "ENTREZID"] <- 8349
genes.df[which(genes.df$SYMBOL == "SDPR"), "ENTREZID"] <- 8436
genes.df[which(genes.df$SYMBOL == "FAM198B"), "ENTREZID"] <- 51313
genes.df[which(genes.df$SYMBOL == "C5orf56"), "ENTREZID"] <- 441108
genes.df[which(genes.df$SYMBOL == "FAM65B"), "ENTREZID"] <- 9750
genes.df[which(genes.df$SYMBOL == "HIST1H1C"), "ENTREZID"] <- 3006
genes.df[which(genes.df$SYMBOL == "HIST1H2AC"), "ENTREZID"] <- 8334
genes.df[which(genes.df$SYMBOL == "HIST1H2BD"), "ENTREZID"] <- 3017
genes.df[which(genes.df$SYMBOL == "C8orf4"), "ENTREZID"] <- 56892
genes.df[which(genes.df$SYMBOL == "FAM214B"), "ENTREZID"] <- 80256
genes.df[which(genes.df$SYMBOL == "C10orf54"), "ENTREZID"] <- 64115
genes.df[which(genes.df$SYMBOL == "OBFC1"), "ENTREZID"] <- 79991
genes.df[which(genes.df$SYMBOL == "RARRES3"), "ENTREZID"] <- 5920
genes.df[which(genes.df$SYMBOL == "LINC00346"), "ENTREZID"] <- 283487
genes.df[which(genes.df$SYMBOL == "FAM214A"), "ENTREZID"] <- 56204
genes.df[which(genes.df$SYMBOL == "FAM65C"), "ENTREZID"] <- 140876

# Remove genes that do not have corresponding Entrez IDs
genes.df <- genes.df[which(!is.na(genes.df$ENTREZID)),]

print(paste0(
  length(genes.df$SYMBOL), "/",
  length(rownames(bj.repro)),
  " RIS markers had Entrez IDs."
))

# Load background gene set
bg.genes <- read.delim("repro_geneIDs_background.tsv")

# GO/KEGG
ris.all.go <- enrichGO(
  gene = as.character(genes.df$ENTREZID),
  universe = as.character(bg.genes$ENTREZID),
  OrgDb = org.Hs.eg.db,
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)

ris.all.kegg <- enrichKEGG(
  gene = as.character(genes.df$ENTREZID),
  universe = as.character(bg.genes$ENTREZID),
  organism = "hsa"
)
ris.all.kegg <- setReadable(
  x = ris.all.kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# Plots of GO/KEGG analysis
dp.go <- dotplot(ris.all.go, showCategory = 20) + ggtitle("GO Biological Processes")
dp.kegg <- dotplot(ris.all.kegg, showCategory = 20) + ggtitle("KEGG Pathways")

png(
  filename = "../../Images/SenotypeComparisons/RIS_GO_all.png",
  width = 1050, height = 1800, res = 160
)
dp.go
dev.off()

png(
  filename = "../../Images/SenotypeComparisons/RIS_KEGG_all.png",
  width = 1000, height = 1400, res = 160
)
dp.kegg
dev.off()


## (RIS markers specific to RIS)
# Convert gene names to Entrez IDs
genes.df.spec <- bitr(
  specific.genes$Gene, 
  fromType = "SYMBOL",
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db, 
  drop = FALSE
)
# Manually input Entrez IDs of genes that failed to map due to alternate names
genes.df.spec[which(genes.df.spec$SYMBOL == "LINC01137"), "ENTREZID"] <- 728431
genes.df.spec[which(genes.df.spec$SYMBOL == "KIAA1324"), "ENTREZID"] <- 57535
genes.df.spec[which(genes.df.spec$SYMBOL == "AHSA2"), "ENTREZID"] <- 130872
genes.df.spec[which(genes.df.spec$SYMBOL == "LINC00152"), "ENTREZID"] <- 112597
genes.df.spec[which(genes.df.spec$DYMBOL == "GPR1"), "ENTREZID"] <- 2825
genes.df.spec[which(genes.df.spec$SYMBOL == "FLJ43879"), "ENTREZID"] <- 401039
genes.df.spec[which(genes.df.spec$SYMBOL == "KIAA1109"), "ENTREZID"] <- 84162
genes.df.spec[which(genes.df.spec$SYMBOL == "SEPT8"), "ENTREZID"] <- 23176
genes.df.spec[which(genes.df.spec$SYMBOL == "FAM26E"), "ENTREZID"] <- 254228
genes.df.spec[which(genes.df.spec$SYMBOL == "FAM188B"), "ENTREZID"] <- 84182
genes.df.spec[which(genes.df.spec$SYMBOL == "ABHD11-AS1"), "ENTREZID"] <- 171022
genes.df.spec[which(genes.df.spec$SYMBOL == "FAM71F2"), "ENTREZID"] <- 346653
genes.df.spec[which(genes.df.spec$SYMBOL == "CMB9-22P13.1"), "ENTREZID"] <- 101927789
genes.df.spec[which(genes.df.spec$SYMBOL == "C11orf70"), "ENTREZID"] <- 85016
genes.df.spec[which(genes.df.spec$SYMBOL == "LYRM5"), "ENTREZID"] <- 144363
genes.df.spec[which(genes.df.spec$SYMBOL == "METTL7A"), "ENTREZID"] <- 25840
genes.df.spec[which(genes.df.spec$SYMBOL == "LHFP"), "ENTREZID"] <- 10186
genes.df.spec[which(genes.df.spec$SYMBOL == "C15orf52"), "ENTREZID"] <- 388115
genes.df.spec[which(genes.df.spec$SYMBOL == "C16orf45"), "ENTREZID"] <- 89927
genes.df.spec[which(genes.df.spec$SYMBOL == "FAM134C"), "ENTREZID"] <- 162427
genes.df.spec[which(genes.df.spec$SYMBOL == "C17orf82"), "ENTREZID"] <- 388407
genes.df.spec[which(genes.df.spec$SYMBOL == "BAIAP2-AS1"), "ENTREZID"] <- 440465
genes.df.spec[which(genes.df.spec$SYMBOL == "ZADH2"), "ENTREZID"] <- 284273
genes.df.spec[which(genes.df.spec$SYMBOL == "MMP24-AS1"), "ENTREZID"] <- 101410538
genes.df.spec[which(genes.df.spec$SYMBOL == "ATP5E"), "ENTREZID"] <- 514
genes.df.spec[which(genes.df.spec$SYMBOL == "MARCH2"), "ENTREZID"] <- 51257
# Remove genes that do not have corresponding Entrez IDs
genes.df.spec <- genes.df.spec[which(!is.na(genes.df.spec$ENTREZID)),]

print(paste0(
  length(genes.df.spec$SYMBOL), "/", 
  length(specific.genes$Gene), 
  " RIS-specific genes had Entrez IDs."
))

# GO/KEGG
ris.spec.go <- enrichGO(
  gene = as.character(genes.df.spec$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)

ris.spec.kegg <- enrichKEGG(
  gene = as.character(genes.df.spec$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  organism = "hsa", 
  qvalueCutoff = 0.05
)
ris.spec.kegg <- setReadable(
  x = ris.spec.kegg, 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID"
)

# Plots of GO/KEGG analysis
dp.go <- dotplot(ris.spec.go, showCategory = 20) + ggtitle("GO Biological Processes")
dp.kegg <- dotplot(ris.spec.kegg, showCategory = 20) + ggtitle("KEGG Pathways")

png(
  filename = "../../Images/SenotypeComparisons/RIS_GO_specific.png", 
  width = 1050, height = 1800, res = 160
)
dp.go
dev.off()

png(
  filename = "../../Images/SenotypeComparisons/RIS_KEGG_specific.png", 
  width = 1050, height = 550, res = 160
)
dp.kegg
dev.off()


## (RIS markers common to other senotypes)
# Convert gene names to Entrez IDs
genes.df.common <- bitr(
  setdiff(rownames(bj.repro), specific.genes$Gene), 
  fromType = "SYMBOL",
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db, 
  drop = FALSE
)
# Manually input Entrez IDs of genes that failed to map due to alternate names
genes.df.common[which(genes.df.common$SYMBOL == "WDR63"), "ENTREZID"] <- 126820
genes.df.common[which(genes.df.common$SYMBOL == "HIST2H2BE"), "ENTREZID"] <- 8349
genes.df.common[which(genes.df.common$SYMBOL == "SDPR"), "ENTREZID"] <- 8436
genes.df.common[which(genes.df.common$SYMBOL == "FAM198B"), "ENTREZID"] <- 51313
genes.df.common[which(genes.df.common$SYMBOL == "C5orf56"), "ENTREZID"] <- 441108
genes.df.common[which(genes.df.common$SYMBOL == "FAM65B"), "ENTREZID"] <- 9750
genes.df.common[which(genes.df.common$SYMBOL == "HIST1H1C"), "ENTREZID"] <- 3006
genes.df.common[which(genes.df.common$SYMBOL == "HIST1H2AC"), "ENTREZID"] <- 8334
genes.df.common[which(genes.df.common$SYMBOL == "HIST1H2BD"), "ENTREZID"] <- 3017
genes.df.common[which(genes.df.common$SYMBOL == "C8orf4"), "ENTREZID"] <- 56892
genes.df.common[which(genes.df.common$SYMBOL == "FAM214B"), "ENTREZID"] <- 80256
genes.df.common[which(genes.df.common$SYMBOL == "C10orf54"), "ENTREZID"] <- 64115
genes.df.common[which(genes.df.common$SYMBOL == "OBFC1"), "ENTREZID"] <- 79991
genes.df.common[which(genes.df.common$SYMBOL == "RARRES3"), "ENTREZID"] <- 5920
genes.df.common[which(genes.df.common$SYMBOL == "LINC00346"), "ENTREZID"] <- 283487
genes.df.common[which(genes.df.common$SYMBOL == "FAM214A"), "ENTREZID"] <- 56204
genes.df.common[which(genes.df.common$SYMBOL == "FAM65C"), "ENTREZID"] <- 140876
# Remove genes that do not have corresponding Entrez IDs
genes.df.common <- genes.df.common[which(!is.na(genes.df.common$ENTREZID)),]

print(paste0(
  length(genes.df.common$SYMBOL), "/", 
  length(setdiff(rownames(bj.repro), specific.genes$Gene)), 
  " genes had Entrez IDs."
))

# GO/KEGG
ris.common.go <- enrichGO(
  gene = as.character(genes.df.common$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)

ris.common.kegg <- enrichKEGG(
  gene = as.character(genes.df.common$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  organism = "hsa", 
  qvalueCutoff = 0.05
)
ris.common.kegg <- setReadable(
  x = ris.common.kegg, 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID"
)

# Plots of GO/KEGG analysis
dp.go <- dotplot(ris.common.go, showCategory = 20) + ggtitle("GO Biological Processes")
dp.kegg <- dotplot(ris.common.kegg, showCategory = 20) + ggtitle("KEGG Pathways")

png(
  filename = "../../Images/SenotypeComparisons/RIS_GO_non-specific.png", 
  width = 1050, height = 1800, res = 160
)
dp.go
dev.off()

png(
  filename = "../../Images/SenotypeComparisons/RIS_KEGG_non-specific.png", 
  width = 1050, height = 650, res = 160
)
dp.kegg
dev.off()



