# Set working directory
setwd("Data/Schiebinger2019/")

# Load packages and Seurat object
library(Seurat)
library(ggplot2)
library(ggrastr)

scbg <- readRDS("scbg_final.rds")

# Identify cells that express high CDKN1A and CDKN2A
scbg$Expression <- "Neither"
exp_cdkn1a <- WhichCells(scbg, expression = Cdkn1a > 2)
exp_cdkn2a <- WhichCells(scbg, expression = Cdkn2a > 1)
exp_both <- intersect(exp_cdkn1a, exp_cdkn2a)

scbg$Expression[which(colnames(scbg) %in% exp_cdkn1a)] <- "CDKN1A"
scbg$Expression[which(colnames(scbg) %in% exp_cdkn2a)] <- "CDKN2A"
scbg$Expression[which(colnames(scbg) %in% exp_both)] <- "Both"

dp.expression <- DimPlot(
  object = scbg, reduction = "umap", 
  group.by = "Expression", raster = TRUE, 
  order = c("Both", "CDKN2A", "CDKN1A", "Neither"), 
  cols = c("#BEBEBE", "#F8766D", "#00BFC4", "#C77CFF")
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Mouse"
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)
dp.expression.split <- DimPlot(
  object = scbg, reduction = "umap", 
  group.by = "Expression", split.by = "Expression", 
  ncol = 2, raster = TRUE, 
  order = c("Both", "CDKN2A", "CDKN1A", "Neither"), 
  cols = c("#BEBEBE", "#F8766D", "#00BFC4", "#C77CFF")
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Mouse"
) + theme(
  text = element_text(size = 36), 
  axis.text = element_blank(), 
  axis.ticks = element_blank()
) + NoLegend()

png(
  filename = "../../Images/MouseReprogramming/scbg_UMAP_expressedgenes.png", 
  width = 1400, height = 1000, res = 160
)
dp.expression
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_UMAP_expressedgenes_split.png", 
  width = 2050, height = 2000, res = 160
)
dp.expression.split
dev.off()

# Markers of CDKN1A-high and CDKN2A-high cells
markers <- FindAllMarkers(
  object = scbg, assay = "RNA", 
  group.by = "Expression", 
  logfc.threshold = 0.5, min.pct = 0.30, 
  only.pos = TRUE
)
markers <- subset(markers, p_val_adj < 0.05)
markers <- markers[order(markers$avg_log2FC, decreasing = TRUE),]

markers.CDKN1A <- subset(markers, cluster == "CDKN1A")
markers.CDKN2A <- subset(markers, cluster == "CDKN2A")
markers.Both <- subset(markers, cluster == "Both")
gc()

print(paste(
  length(unique(markers.CDKN1A$gene)), 
  "genes were identified as markers for CDKN1A-high cells."
))
print(paste(
  length(unique(markers.CDKN2A$gene)), 
  "genes were identified as markers for CDKN2A-high cells."
))
print(paste(
  length(unique(markers.Both$gene)), 
  "genes were identified as markers for double-high cells."
))

# Search for SenMayo genes within markers
senmayo <- read.csv("../../Tables/GeneSets/senmayo_genes_mouse.csv")
sm.Both <- markers.Both[which(markers.Both$gene %in% senmayo$Gene.murine.),]
sm.CDKN1A <- markers.CDKN1A[which(markers.CDKN1A$gene %in% senmayo$Gene.murine.),]
sm.CDKN2A <- markers.CDKN2A[which(markers.CDKN2A$gene %in% senmayo$Gene.murine.),]

print(paste0(
  length(unique(c(sm.CDKN1A$gene, sm.CDKN2A$gene, sm.Both$gene))), "/", 
  length(unique(c(markers.CDKN1A$gene, markers.CDKN2A$gene, markers.Both$gene))), 
  " marker genes were identified as being part of the SenMayo gene set."
))
print(paste(
  length(unique(sm.CDKN1A$gene)), 
  "SenMayo genes were markers of CDKN1A-high cells."
))
print(paste(
  length(unique(sm.CDKN2A$gene)), 
  "SenMayo genes were markers of CDKN2A-high cells."
))
print(paste(
  length(unique(sm.Both$gene)), 
  "SenMayo genes were markers of double-high cells."
))

print(paste(
  "There are", 
  length(intersect(unique(sm.Both$gene), unique(sm.CDKN1A$gene))), 
  "SenMayo genes common to both CDKN1A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(sm.Both$gene), unique(sm.CDKN2A$gene))), 
  "SenMayo genes common to both CDKN2A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(sm.CDKN1A$gene), unique(sm.CDKN2A$gene))), 
  "SenMayo genes common to both CDKN1A-high and CDKN2A-high cells."
))
print(paste(
  "There are", 
  length(Reduce(intersect, list(unique(sm.Both$gene), unique(sm.CDKN1A$gene), unique(sm.CDKN2A$gene)))), 
  "SenMayo genes common to all three groups."
))

# Create table of SenMayo markers in dataset
senmayo.markers <- data.frame(matrix(
  ncol = 4, 
  nrow = length(unique(c(sm.Both$gene, sm.CDKN1A$gene, sm.CDKN2A$gene)))
))
colnames(senmayo.markers) <- c("Gene", "CDKN1A-high", "CDKN2A-high", "Double-high")
senmayo.markers$Gene <- unique(c(sm.Both$gene, sm.CDKN1A$gene, sm.CDKN2A$gene))

for (gene in senmayo.markers$Gene) {
  senmayo.markers[which(senmayo.markers$Gene == gene), "CDKN1A-high"] <- ifelse(gene %in% sm.CDKN1A$gene, "Present", "Absent")
  senmayo.markers[which(senmayo.markers$Gene == gene), "CDKN2A-high"] <- ifelse(gene %in% sm.CDKN2A$gene, "Present", "Absent")
  senmayo.markers[which(senmayo.markers$Gene == gene), "Double-high"] <- ifelse(gene %in% sm.Both$gene, "Present", "Absent")
}

write.table(
  senmayo.markers, file = "../../Tables/MouseReprogramming/senmayo_markers_schiebinger.tsv", 
  sep = "\t", row.names = FALSE, qmethod = "double"
)


## SenSig
# Create list of relevant genes from SenSig
ss.filepaths = c(
  "../../Tables/MouseReprogramming/GeneSets/SenSig_ligands_mouse_VMLscRNAseqFibrotic.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Common-P1.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Common-P2.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Common-P3.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Fibrotic-P.csv",
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Fibrotic.csv",
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Myo1-P.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Myo1.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_PeriMyo2-P.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Myo2-P.csv", 
  "../../Tables/MouseReprogramming/GeneSets/SenSig_genes_sortedcells_Myo2.csv"
)
gene.names <- vector(mode = "list", length = length(ss.filepaths))

for (i in 1:length(ss.filepaths)) {
  # Common-P2 file has duplicate row names
  if (grepl("Common-P2", ss.filepaths[i])) {
    row.vector <- NULL
  } else {
    row.vector <- 1
  }
  
  # Download SenSig gene table
  senfile <- read.csv(file = ss.filepaths[i], row.names = row.vector)
  senfile <- subset(senfile, p_val_adj < 0.05)
  
  # Obtain gene names from table
  if ("gene" %in% colnames(senfile)) {
    genes <- senfile$gene
  } else {
    genes <- rownames(senfile)
  }
  
  gene.names[[i]] <- genes
}

sensig.genes <- unique(unlist(gene.names))

# Search for SenSig genes within markers
ss.Both <- markers.Both[which(markers.Both$gene %in% sensig.genes),]
ss.CDKN1A <- markers.CDKN1A[which(markers.CDKN1A$gene %in% sensig.genes),]
ss.CDKN2A <- markers.CDKN2A[which(markers.CDKN2A$gene %in% sensig.genes),]

print(paste0(
  length(unique(c(ss.CDKN1A$gene, ss.CDKN2A$gene, ss.Both$gene))), "/", 
  length(unique(c(markers.CDKN1A$gene, markers.CDKN2A$gene, markers.Both$gene))), 
  " marker genes were identified as being part of the SenSig gene set."
))
print(paste(
  length(ss.CDKN1A$gene), 
  "SenSig genes were markers of CDKN1A-high cells."
))
print(paste(
  length(ss.CDKN2A$gene), 
  "SenSig genes were markers of CDKN2A-high cells."
))
print(paste(
  length(ss.Both$gene), 
  "SenSig genes were markers of double-high cells."
))

print(paste(
  "There are", 
  length(intersect(unique(ss.Both$gene), unique(ss.CDKN1A$gene))), 
  "SenSig genes common to both CDKN1A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(ss.Both$gene), unique(ss.CDKN2A$gene))), 
  "SenSig genes common to both CDKN2A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(ss.CDKN1A$gene), unique(ss.CDKN2A$gene))), 
  "SenSig genes common to both CDKN1A-high and CDKN2A-high cells."
))
print(paste(
  "There are", 
  length(Reduce(intersect, list(unique(ss.Both$gene), unique(ss.CDKN1A$gene), unique(ss.CDKN2A$gene)))), 
  "SenSig genes common to all three groups."
))

# Create table of SenSig markers in dataset
sensig.markers <- data.frame(matrix(
  ncol = 4, 
  nrow = length(unique(c(ss.Both$gene, ss.CDKN1A$gene, ss.CDKN2A$gene)))
))
colnames(sensig.markers) <- c("Gene", "CDKN1A-high", "CDKN2A-high", "Double-high")
sensig.markers$Gene <- unique(c(ss.Both$gene, ss.CDKN1A$gene, ss.CDKN2A$gene))

for (gene in sensig.markers$Gene) {
  sensig.markers[which(sensig.markers$Gene == gene), "CDKN1A-high"] <- ifelse(gene %in% ss.CDKN1A$gene, "Present", "Absent")
  sensig.markers[which(sensig.markers$Gene == gene), "CDKN2A-high"] <- ifelse(gene %in% ss.CDKN2A$gene, "Present", "Absent")
  sensig.markers[which(sensig.markers$Gene == gene), "Double-high"] <- ifelse(gene %in% ss.Both$gene, "Present", "Absent")
}

write.table(
  sensig.markers, file = "../../Tables/MouseReprogramming/sensig_markers_schiebinger.tsv", 
  sep = "\t", row.names = FALSE, qmethod = "double"
)



