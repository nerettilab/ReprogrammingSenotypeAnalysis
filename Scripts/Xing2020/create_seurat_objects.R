# Set working directory
setwd("Data/Xing2020")

# Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)

# Load RNA data
counts <- read.csv(
  file = "GSE118258/GSE118258_UMI.csv", 
  sep = ",", row.names = 1
)

# Obtain gene names from Ensembl IDs
idmap <- select(
  x = EnsDb.Hsapiens.v86, 
  filter = GeneIdFilter("ENS", "startsWith"), 
  keys = rownames(counts), 
  column = "SYMBOL", 
  keytype = "GENEID"
)

# Modify gene symbols to get rid of duplicate names
idmap <- idmap[-c(which(idmap$GENEID == "ENSG00000271858"), 
         which(idmap$GENEID == "ENSG00000235271"), 
         which(idmap$GENEID == "ENSG00000257613"), 
         which(idmap$GENEID == "ENSG00000234229"), 
         which(idmap$GENEID == "ENSG00000231963"), 
         which(idmap$GENEID == "ENSG00000225655"), 
         which(idmap$GENEID == "ENSG00000272167"), 
         which(idmap$GENEID == "ENSG00000232995"), 
         which(idmap$GENEID == "ENSG00000273259")),]
idmap[which(idmap$GENEID == "ENSG00000272617"),]$SYMBOL <- "COG8-PDF"
idmap[which(idmap$GENEID == "ENSG00000227540"),]$SYMBOL <- "RP11-152N13.5"
idmap[which(idmap$GENEID == "ENSG00000229694"),]$SYMBOL <- "RP11-305L7.6"
idmap[which(idmap$GENEID == "ENSG00000231512"),]$SYMBOL <- "SEPTIN14P21"
idmap[which(idmap$GENEID == "ENSG00000257815"),]$SYMBOL <- "PRANCR"
idmap[which(idmap$GENEID == "ENSG00000268154"),]$SYMBOL <- "Metazoa_SRP.25"
idmap[which(idmap$GENEID == "ENSG00000203286"),]$SYMBOL <- "Metazoa_SRP.3"
idmap[which(idmap$GENEID == "ENSG00000214305"),]$SYMBOL <- "Metazoa_SRP.8"
idmap[which(idmap$GENEID == "ENSG00000269103"),]$SYMBOL <- "Metazoa_SRP.26"
idmap[which(idmap$GENEID == "ENSG00000241111"),]$SYMBOL <- "RP11-129B22.1"
idmap[which(idmap$GENEID == "ENSG00000268592"),]$SYMBOL <- "RP11-244K5.8"
idmap[which(idmap$GENEID == "ENSG00000228741"),]$SYMBOL <- "SPATA13-LncRNA"
idmap[which(idmap$GENEID == "ENSG00000255104"),]$SYMBOL <- "ZNF286A-TBC1D26"
idmap[which(idmap$GENEID == "ENSG00000187838"),]$SYMBOL <- "PLSCR3"

# Replace Ensembl IDs with gene names when possible
rownames(counts) <- replace(rownames(counts), rownames(counts) %in% idmap$GENEID, idmap$SYMBOL)

# Create Seurat object
xing <- CreateSeuratObject(
  counts = counts,  
  assay = "RNA", 
  project = "Xing", 
  min.cells = 3, min.features = 200
)
xing$dataset <- "Xing"

# Annotate cells with their sample of origin
annotations <- read.table(
  file = "GSE118258/GSE118258_Annotation.txt", 
  col.names = c("Cell_ID", "Timepoint"), 
  skip = 1
)
xing@meta.data$Timepoint <- annotations$Timepoint

Idents(xing) <- "Timepoint"
xing <- RenameIdents(
  object = xing, 
  "D0" = "D00", 
  "D2" = "D02", 
  "D8" = "D08", 
  "D12" = "D12", 
  "D16_negative" = "D16", 
  "D16_positive" = "D16"
)
xing[["Timepoint"]] <- Idents(xing)

xing[["percent.mt"]] <- PercentageFeatureSet(
  object = xing, 
  pattern = "^MT-"
)

saveRDS(xing, file = "xing.rds")

# Garbage collection
rm(counts, idmap, annotations)
gc(reset=TRUE)

# Quality Control
VlnPlot(
  object = xing, 
  features = c("nFeature_RNA", "nCount_RNA"), 
  pt.size = 0.1, 
  group.by = "Timepoint", 
  ncol = 2, 
  layer = "counts"
)
FeatureScatter(
  object = xing, 
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
  group.by = "Timepoint"
) + NoLegend()
FeatureScatter(
  object = xing, 
  feature1 = "nCount_RNA", feature2 = "percent.mt", 
  group.by = "Timepoint"
) + geom_hline(yintercept = 10)

# Filter out low-quality cells
xing.list <- SplitObject(xing, split.by = "Timepoint")

for (i in 1:length(xing.list)) {
  xing.list[[i]] <- subset(
    x = xing.list[[i]], 
    subset = nFeature_RNA > quantile(xing.list[[i]]$nFeature_RNA, probs = 0.05) & nFeature_RNA < quantile(xing.list[[i]]$nFeature_RNA, probs = 0.95) & percent.mt < 10
  )
}

xing <- merge(
  x = xing.list[[1]], 
  y = xing.list[-1], 
  project = "Xing"
)
xing <- JoinLayers(xing, assay = "RNA")
Idents(xing) <- "Timepoint"

saveRDS(xing, file = "xing_filtered.rds")

rm(xing.list, i)
gc(reset = TRUE)

# Process merged Seurat object
xing <- NormalizeData(xing)
xing <- FindVariableFeatures(xing)
xing <- ScaleData(xing, vars.to.regress = "percent.mt")
xing <- RunPCA(xing)
xing <- RunUMAP(xing, dims = 1:30)
xing <- FindNeighbors(xing, dims = 1:30)
xing <- FindClusters(xing, resolution = 0.5)
DimPlot(xing, reduction = "umap", group.by = "Timepoint", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP by Treatment Time")

saveRDS(xing, file = "xing_final.rds")


