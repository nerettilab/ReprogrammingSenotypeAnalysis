# Set working directory
setwd("Data/MergedHuman/")

# Load packages and Seurat object
library(Seurat)
library(ggplot2)
library(biomaRt)

combined <- readRDS("combined.rds")


# Identify cells that express high CDKN1A and CDKN2A
combined$Expression <- "Neither"
exp_cdkn1a <- WhichCells(combined, expression = CDKN1A > 1)
exp_cdkn2a <- WhichCells(combined, expression = CDKN2A > 1)
exp_both <- intersect(exp_cdkn1a, exp_cdkn2a)

combined$Expression[which(colnames(combined) %in% exp_cdkn1a)] <- "CDKN1A"
combined$Expression[which(colnames(combined) %in% exp_cdkn2a)] <- "CDKN2A"
combined$Expression[which(colnames(combined) %in% exp_both)] <- "Both"

dp.expression <- DimPlot(
  object = combined, reduction = "umap", 
  group.by = "Expression", 
  order = c("Both", "CDKN2A", "CDKN1A", "Neither"), 
  cols = c("#BEBEBE", "#F8766D", "#00BFC4", "#C77CFF")
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Human"
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)
dp.expression.split <- DimPlot(
  object = combined, reduction = "umap", 
  group.by = "Expression", split.by = "Expression", 
  order = c("Both", "CDKN2A", "CDKN1A", "Neither"), 
  cols = c("#BEBEBE", "#F8766D", "#00BFC4", "#C77CFF"), 
  ncol = 2
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Human"
) + theme(
  text = element_text(size = 36), 
  axis.text = element_blank(), 
  axis.ticks = element_blank()
) + NoLegend()

png(
  filename = "../../Images/HumanReprogramming/xing+pana_UMAP_expressedgenes.png", 
  width = 1400, height = 1000, res = 160
)
dp.expression
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_UMAP_expressedgenes_split.png", 
  width = 2020, height = 1780, res = 160
)
dp.expression.split
dev.off()

# Markers of CDKN1A-high, CDKN2A-high, and double-high cells
Idents(combined) <- "Expression"
for (category in c("CDKN1A", "CDKN2A", "Both")) {
  # Only consider genes found in ≥30% of cells in category of interest
  # Only consider genes with a log2FC ≥ 0.5
  markers <- FindConservedMarkers(
    combined, ident.1 = category, 
    grouping.var = "dataset", only.pos = TRUE, 
    logfc.threshold = 0.5, min.pct = 0.30
  )
  
  markers <- subset(markers, Xing_p_val_adj < 0.05 & Panariello_p_val_adj < 0.05)
  markers$Xing_p_val <- NULL
  markers$Panariello_p_val <- NULL
  markers$max_pval <- NULL
  markers$minimump_p_val <- NULL
  
  # Order cells by the average of the log2FC values in both datasets
  markers$AVG_avg_log2FC <- (markers$Xing_avg_log2FC + markers$Panariello_avg_log2FC) / 2
  markers <- markers[order(markers$AVG_avg_log2FC, decreasing = TRUE),]
  
  assign(paste0("markers.", category), markers)
  gc()
}


print(paste(
  length(rownames(markers.CDKN1A)), 
  "genes were identified as markers for CDKN1A-high cells."
))
print(paste(
  length(rownames(markers.CDKN2A)), 
  "genes were identified as markers for CDKN2A-high cells."
))
print(paste(
  length(rownames(markers.Both)), 
  "genes were identified as markers for double-high cells."
))
print(paste(
  "There were a total of", 
  length(unique(c(rownames(markers.Both), rownames(markers.CDKN1A), rownames(markers.CDKN2A)))), 
  "unique identified markers."
))


## SenMayo
# Search for SenMayo genes within markers
senmayo <- read.csv("../../Tables/HumanReprogramming/GeneSets/senmayo_genes_human.csv")
sm.Both <- markers.Both[which(rownames(markers.Both) %in% senmayo$Gene.human.),]
sm.CDKN1A <- markers.CDKN1A[which(rownames(markers.CDKN1A) %in% senmayo$Gene.human.),]
sm.CDKN2A <- markers.CDKN2A[which(rownames(markers.CDKN2A) %in% senmayo$Gene.human.),]

print(paste0(
  length(unique(c(rownames(sm.CDKN1A), rownames(sm.CDKN2A), rownames(sm.Both)))), "/", 
  length(unique(c(rownames(markers.CDKN1A), rownames(markers.CDKN2A), rownames(markers.Both)))), 
  " marker genes were identified as being part of the SenMayo gene set."
))
print(paste(
  length(rownames(sm.CDKN1A)), 
  "SenMayo genes were markers of CDKN1A-high cells."
))
print(paste(
  length(rownames(sm.CDKN2A)), 
  "SenMayo genes were markers of CDKN2A-high cells."
))
print(paste(
  length(rownames(sm.Both)), 
  "SenMayo genes were markers of double-high cells."
))

print(paste(
  "There are", 
  length(intersect(unique(rownames(sm.Both)), unique(rownames(sm.CDKN1A)))), 
  "SenMayo genes common to both CDKN1A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(rownames(sm.Both)), unique(rownames(sm.CDKN2A)))), 
  "SenMayo genes common to both CDKN2A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(rownames(sm.CDKN1A)), unique(rownames(sm.CDKN2A)))), 
  "SenMayo genes common to both CDKN1A-high and CDKN2A-high cells."
))
print(paste(
  "There are", 
  length(Reduce(intersect, list(unique(rownames(sm.Both)), unique(rownames(sm.CDKN1A)), unique(rownames(sm.CDKN2A))))), 
  "SenMayo genes common to all three groups."
))

# Create table of SenMayo markers in dataset
senmayo.markers <- data.frame(matrix(
  ncol = 4, 
  nrow = length(unique(c(rownames(sm.Both), rownames(sm.CDKN1A), rownames(sm.CDKN2A))))
))
colnames(senmayo.markers) <- c("Gene", "CDKN1A-high", "CDKN2A-high", "Double-high")
senmayo.markers$Gene <- unique(c(rownames(sm.Both), rownames(sm.CDKN1A), rownames(sm.CDKN2A)))

for (gene in senmayo.markers$Gene) {
  senmayo.markers[which(senmayo.markers$Gene == gene), "CDKN1A-high"] <- ifelse(gene %in% rownames(sm.CDKN1A), "Present", "Absent")
  senmayo.markers[which(senmayo.markers$Gene == gene), "CDKN2A-high"] <- ifelse(gene %in% rownames(sm.CDKN2A), "Present", "Absent")
  senmayo.markers[which(senmayo.markers$Gene == gene), "Double-high"] <- ifelse(gene %in% rownames(sm.Both), "Present", "Absent")
}

write.table(
  senmayo.markers, file = "../../Tables/HumanReprogramming/senmayo_markers_xing+pana.tsv", 
  sep = "\t", row.names = FALSE, qmethod = "double"
)


## SenSig
# Download human and mouse biomaRt databases
human.mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse.mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Create list of relevant genes from SenSig
ss.filepaths = c(
  "../../Tables/HumanReprogramming/GeneSets/SenSig_ligands_human_IPFscRNAseqKRT5-KRT17+.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_ligands_human_IPFscRNAseqMyofibroblast.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_ligands_human_CancerscRNAseqCAFs.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_ligands_human_CapsulescRNAseqCluster2.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_ligands_human_CapsulescRNAseqCluster3.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Common-P1.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Common-P2.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Common-P3.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Fibrotic-P.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Fibrotic.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Myo1-P.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Myo1.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_PeriMyo2-P.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Myo2-P.csv", 
  "../../Tables/HumanReprogramming/GeneSets/SenSig_genes_sortedcells_Myo2.csv"
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
  
  ## Obtain gene names from table
  if ("gene" %in% colnames(senfile)) {
    
    genes <- senfile$gene
    
    # Obtain Ensembl IDs of mouse genes
    bm.mouse <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "mgi_symbol"), 
      filters = "external_gene_name", 
      values = genes, 
      mart = mouse.mart
    )
    
    # Retrieve Ensembl IDs of homologous human genes
    homologs <- getHomologs(bm.mouse$ensembl_gene_id, "mus_musculus", "homo_sapiens")
    homolog.ids <- unique(homologs$hsapiens_homolog_ensembl_gene)
    homolog.ids <- homolog.ids[which(homolog.ids != "")]
    
    # Convert Ensembl IDs to human gene symbols
    bm.human <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"), 
      filters = "ensembl_gene_id", 
      values = homolog.ids, 
      mart = human.mart
    )
    genes <- unique(bm.human$external_gene_name)
    genes <- genes[which(genes != "")]
    
  } else {
    
    genes <- rownames(senfile)
    
  }
  
  gene.names[[i]] <- genes
}

sensig.genes <- unique(unlist(gene.names))

# Search for SenSig genes within markers
ss.Both <- markers.Both[which(rownames(markers.Both) %in% sensig.genes),]
ss.CDKN1A <- markers.CDKN1A[which(rownames(markers.CDKN1A) %in% sensig.genes),]
ss.CDKN2A <- markers.CDKN2A[which(rownames(markers.CDKN2A) %in% sensig.genes),]

print(paste0(
  length(unique(c(rownames(ss.CDKN1A), rownames(ss.CDKN2A), rownames(ss.Both)))), "/", 
  length(unique(c(rownames(markers.CDKN1A), rownames(markers.CDKN2A), rownames(markers.Both)))), 
  " marker genes were identified as being part of the SenSig gene set."
))
print(paste(
  length(rownames(ss.CDKN1A)), 
  "SenSig genes were markers of CDKN1A-high cells."
))
print(paste(
  length(rownames(ss.CDKN2A)), 
  "SenSig genes were markers of CDKN2A-high cells."
))
print(paste(
  length(rownames(ss.Both)), 
  "SenSig genes were markers of double-high cells."
))

print(paste(
  "There are", 
  length(intersect(unique(rownames(ss.Both)), unique(rownames(ss.CDKN1A)))), 
  "SenSig genes common to both CDKN1A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(rownames(ss.Both)), unique(rownames(ss.CDKN2A)))), 
  "SenSig genes common to both CDKN2A-high and double-high cells."
))
print(paste(
  "There are", 
  length(intersect(unique(rownames(ss.CDKN1A)), unique(rownames(ss.CDKN2A)))), 
  "SenSig genes common to both CDKN1A-high and CDKN2A-high cells."
))
print(paste(
  "There are", 
  length(Reduce(intersect, list(unique(rownames(ss.Both)), unique(rownames(ss.CDKN1A)), unique(rownames(ss.CDKN2A))))), 
  "SenSig genes common to all three groups."
))

# Create table of SenSig markers in dataset
sensig.markers <- data.frame(matrix(
  ncol = 4, 
  nrow = length(unique(c(rownames(ss.Both), rownames(ss.CDKN1A), rownames(ss.CDKN2A))))
))
colnames(sensig.markers) <- c("Gene", "CDKN1A-high", "CDKN2A-high", "Double-high")
sensig.markers$Gene <- unique(c(rownames(ss.Both), rownames(ss.CDKN1A), rownames(ss.CDKN2A)))

for (gene in sensig.markers$Gene) {
  sensig.markers[which(sensig.markers$Gene == gene), "CDKN1A-high"] <- ifelse(gene %in% rownames(ss.CDKN1A), "Present", "Absent")
  sensig.markers[which(sensig.markers$Gene == gene), "CDKN2A-high"] <- ifelse(gene %in% rownames(ss.CDKN2A), "Present", "Absent")
  sensig.markers[which(sensig.markers$Gene == gene), "Double-high"] <- ifelse(gene %in% rownames(ss.Both), "Present", "Absent")
}

write.table(
  sensig.markers, file = "../../Tables/sensig_markers_xing+pana.tsv", 
  sep = "\t", row.names = FALSE, qmethod = "double"
)


