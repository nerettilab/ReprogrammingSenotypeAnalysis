# Set working directory
setwd("Data/Panariello2023/")

# Load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)

## Load data for each timepoint
# Day 0
exp.matrix.D00A <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D00/GSM6894025_D0A_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D00/GSM6894025_D0A_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D00/GSM6894025_D0A_genes.tsv"
)

exp.matrix.D00B <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D00/GSM6894026_D0B_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D00/GSM6894026_D0B_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D00/GSM6894026_D0B_genes.tsv"
)

# Day 3
exp.matrix.D03A <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D03/GSM6894027_D3A_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D03/GSM6894027_D3A_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D03/GSM6894027_D3A_genes.tsv"
)

exp.matrix.D03B <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D03/GSM6894028_D3B_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D03/GSM6894028_D3B_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D03/GSM6894028_D3B_genes.tsv"
)

# Day 5
exp.matrix.D05A1 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894029_D5A1_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894029_D5A1_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894029_D5A1_genes.tsv"
)

exp.matrix.D05A2 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894030_D5A2_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894030_D5A2_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894030_D5A2_genes.tsv"
)

exp.matrix.D05B1 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894031_D5B1_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894031_D5B1_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894031_D5B1_genes.tsv"
)

exp.matrix.D05B2 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894032_D5B2_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894032_D5B2_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D05/GSM6894032_D5B2_genes.tsv"
)

# Day 7
exp.matrix.D07A <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894033_D7A_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894033_D7A_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894033_D7A_genes.tsv"
)

exp.matrix.D07B1 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894034_D7B1_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894034_D7B1_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894034_D7B1_genes.tsv"
)

exp.matrix.D07B2 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894035_D7B2_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894035_D7B2_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D07/GSM6894035_D7B2_genes.tsv"
)

# Day 9
exp.matrix.D09A1 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894036_D9A1_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894036_D9A1_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894036_D9A1_genes.tsv"
)

exp.matrix.D09A2 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894037_D9A2_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894037_D9A2_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894037_D9A2_genes.tsv"
)

exp.matrix.D09B1 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894038_D9B1_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894038_D9B1_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894038_D9B1_genes.tsv"
)

exp.matrix.D09B2 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894039_D9B2_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894039_D9B2_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D09/GSM6894039_D9B2_genes.tsv"
)

# Day 11
exp.matrix.D11A <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894040_D11A_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894040_D11A_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894040_D11A_genes.tsv"
)

exp.matrix.D11B1 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894041_D11B1_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894041_D11B1_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894041_D11B1_genes.tsv"
)

exp.matrix.D11B2 <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894042_D11B2_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894042_D11B2_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D11/GSM6894042_D11B2_genes.tsv"
)

# Day 13
exp.matrix.D13A <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D13/GSM6894043_D13A_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D13/GSM6894043_D13A_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D13/GSM6894043_D13A_genes.tsv"
)

exp.matrix.D13B <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D13/GSM6894044_D13B_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D13/GSM6894044_D13B_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D13/GSM6894044_D13B_genes.tsv"
)

# Day 15
exp.matrix.D15A <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D15/GSM6894045_D15A_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D15/GSM6894045_D15A_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D15/GSM6894045_D15A_genes.tsv"
)

exp.matrix.D15B <- ReadMtx(
  mtx = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D15/GSM6894046_D15B_matrix.mtx", 
  cells = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D15/GSM6894046_D15B_barcodes.tsv", 
  features = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/GSE221739/D15/GSM6894046_D15B_genes.tsv"
)


# Sample list
samples <- c("exp.matrix.D00A", "exp.matrix.D00B", "exp.matrix.D03A", "exp.matrix.D03B", "exp.matrix.D05A1", "exp.matrix.D05A2", "exp.matrix.D05B1", "exp.matrix.D05B2", "exp.matrix.D07A", "exp.matrix.D07B1", "exp.matrix.D07B2", "exp.matrix.D09A1", "exp.matrix.D09A2", "exp.matrix.D09B1", "exp.matrix.D09B2", "exp.matrix.D11A", "exp.matrix.D11B1", "exp.matrix.D11B2", "exp.matrix.D13A", "exp.matrix.D13B", "exp.matrix.D15A", "exp.matrix.D15B")

output_dir <- "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/seurat_objects_raw"

# Create and assign each Seurat object by name
for (s in samples) {
  newname <- str_remove_all(s, "exp.matrix.")
  obj <- get(s)
  
  obj <- CreateSeuratObject(
    counts = obj, project = "Panariello", 
    min.cells = 3, min.features = 200
  )
  obj$orig.ident <- paste0("pana_", newname)
  obj$dataset <- "Panariello"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  assign(newname, obj)  # Save to memory with name like "D00A"
  
  saveRDS(obj, file = file.path(output_dir, paste0(newname, "_raw.rds")))
}

# Garbage collection
rm(exp.matrix.D00A, exp.matrix.D00B, exp.matrix.D03A, 
   exp.matrix.D03B, exp.matrix.D05A1, exp.matrix.D05A2, 
   exp.matrix.D05B1, exp.matrix.D05B2, exp.matrix.D07A, 
   exp.matrix.D07B1, exp.matrix.D07B2, exp.matrix.D09A1, 
   exp.matrix.D09A2, exp.matrix.D09B1, exp.matrix.D09B2, 
   exp.matrix.D11A, exp.matrix.D11B1, exp.matrix.D11B2, 
   exp.matrix.D13A, exp.matrix.D13B, exp.matrix.D15A, 
   exp.matrix.D15B, obj, newname, output_dir)
gc(reset = TRUE)

## Quality Control
# Sample list
samples <- c("D00A", "D00B", "D03A", "D03B", "D05A1", 
             "D05A2", "D05B1", "D05B2", "D07A", "D07B1", 
             "D07B2", "D09A1", "D09A2", "D09B1", "D09B2", 
             "D11A", "D11B1", "D11B2", "D13A", "D13B", 
             "D15A", "D15B")

# Violin Plots
for (s in samples) {
  obj <- get(s) # Get raw object from memory
  print(VlnPlot(
    object = obj, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    layer = "counts"))
}

# Counts vs Features
for (s in samples) {
  obj <- get(s) # Get raw object from memory
  print(FeatureScatter(
    object = obj, 
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
    ) + 
    geom_hline(yintercept = quantile(obj$nFeature_RNA, probs = 0.05)) + 
    geom_hline(yintercept = quantile(obj$nFeature_RNA, probs = 0.95))
  )
}

# Counts vs Mitochondrial Content
for (s in samples) {
  obj <- get(s) # Get raw object from memory
  print(FeatureScatter(
    object = obj, 
    feature1 = "nCount_RNA", feature2 = "percent.mt"
  ) + geom_hline(yintercept = 10))
}

# Filter out low-quality cells
filtered_dir <- "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/seurat_objects_filtered/"

for (s in samples) {
  obj <- get(s)  # Get raw object from memory
  
  # Apply filtering
  obj_filtered <- subset(
    obj, 
    subset = nFeature_RNA > quantile(obj$nFeature_RNA, probs = 0.05) & nFeature_RNA < quantile(obj$nFeature_RNA, probs = 0.95) & percent.mt < 10
  )
  
  # Save filtered version in memory
  assign(s, obj_filtered)
  
  # Save filtered object to disk
  saveRDS(obj_filtered, file = file.path(filtered_dir, paste0(s, "_filtered.rds")))
}

# Merge filtered Seurat objects
seurat_list <- list(D00A, D00B, D03A, D03B, D05A1, D05A2, 
                    D05B1, D05B2, D07A, D07B1, D07B2, D09A1, 
                    D09A2, D09B1, D09B2, D11A, D11B1, D11B2, 
                    D13A, D13B, D15A, D15B)

combined <- merge(
  x = seurat_list[[1]], y = seurat_list[-1], 
  add.cell.ids = samples, 
  project = "Panariello"
)
combined <- JoinLayers(object = combined, assay = "RNA")
Idents(combined) <- "orig.ident"

# Garbage collection
rm(D00A, D00B, D03A, D03B, D05A1, D05A2, D05B1, D05B2, D07A, 
   D07B1, D07B2, D09A1, D09A2, D09B1, D09B2, D11A, D11B1, 
   D11B2, D13A, D13B, D15A, D15B, obj, s, samples, 
   seurat_list)
gc(reset = TRUE)

# Process merged Seurat object
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined, vars.to.regress = "percent.mt")
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
DimPlot(
  object = combined, reduction = "umap", 
  group.by = "orig.ident", label = TRUE, repel = TRUE
) + ggtitle("UMAP by Sample")

# Add time point metadata
metadata <- combined@meta.data
metadata$Timepoint <- NA

metadata[which(metadata$orig.ident == "pana_D00A"),]$Timepoint <- "D00"
metadata[which(metadata$orig.ident == "pana_D00B"),]$Timepoint <- "D00"
metadata[which(metadata$orig.ident == "pana_D03A"),]$Timepoint <- "D03"
metadata[which(metadata$orig.ident == "pana_D03B"),]$Timepoint <- "D03"
metadata[which(metadata$orig.ident == "pana_D05A1"),]$Timepoint <- "D05"
metadata[which(metadata$orig.ident == "pana_D05A2"),]$Timepoint <- "D05"
metadata[which(metadata$orig.ident == "pana_D05B1"),]$Timepoint <- "D05"
metadata[which(metadata$orig.ident == "pana_D05B2"),]$Timepoint <- "D05"
metadata[which(metadata$orig.ident == "pana_D07A"),]$Timepoint <- "D07"
metadata[which(metadata$orig.ident == "pana_D07B1"),]$Timepoint <- "D07"
metadata[which(metadata$orig.ident == "pana_D07B2"),]$Timepoint <- "D07"
metadata[which(metadata$orig.ident == "pana_D09A1"),]$Timepoint <- "D09"
metadata[which(metadata$orig.ident == "pana_D09A2"),]$Timepoint <- "D09"
metadata[which(metadata$orig.ident == "pana_D09B1"),]$Timepoint <- "D09"
metadata[which(metadata$orig.ident == "pana_D09B2"),]$Timepoint <- "D09"
metadata[which(metadata$orig.ident == "pana_D11A"),]$Timepoint <- "D11"
metadata[which(metadata$orig.ident == "pana_D11B1"),]$Timepoint <- "D11"
metadata[which(metadata$orig.ident == "pana_D11B2"),]$Timepoint <- "D11"
metadata[which(metadata$orig.ident == "pana_D13A"),]$Timepoint <- "D13"
metadata[which(metadata$orig.ident == "pana_D13B"),]$Timepoint <- "D13"
metadata[which(metadata$orig.ident == "pana_D15A"),]$Timepoint <- "D15"
metadata[which(metadata$orig.ident == "pana_D15B"),]$Timepoint <- "D15"

combined@meta.data <- metadata

DimPlot(
  object = combined, reduction = "umap", 
  group.by = "Timepoint", label = TRUE, repel = TRUE
) + ggtitle("UMAP by Treatment Time")

# Save final merged Seurat object
saveRDS(combined, file = "/oscar/data/nneretti/tnunez2/SenescenceProgramOSKM2026/Panariello2023/pana_merged_final.rds")


