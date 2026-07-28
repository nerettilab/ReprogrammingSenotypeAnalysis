# Set working directory
setwd("Data/MergedHuman/")

# Load packages and Seurat objects
library(Seurat)
library(ggplot2)

xing <- readRDS("../Xing2020/xing_final.rds")
pana <- readRDS("../Panariello2023/pana_merged_final.rds")

combined <- merge(
  x = xing, y = pana, 
  project = "repro_human"
)

rm(xing, pana)
gc()

# Process merged Seurat object
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined, vars.to.regress = c("percent.mt"))
combined <- RunPCA(combined)

# Integrate based on dataset of origin
options(future.globals.maxSize = 25 * 1000 * 1024^2)
combined <- IntegrateLayers(
  object = combined, 
  method = RPCAIntegration, 
  orig = "pca", new.reduction = "rpca", 
  verbose = FALSE
)
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

gc()

# Continue processing
combined <- FindNeighbors(combined, dims = 1:30, reduction = "rpca")
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:30, reduction = "rpca")

# Annotate time points with the datasets they originate from
Idents(combined) <- "Timepoint"
combined <- RenameIdents(
  object = combined, 
  "D00" = "D00", 
  "D02" = "D02-X", 
  "D03" = "D03-P", 
  "D05" = "D05-P", 
  "D07" = "D07-P", 
  "D08" = "D08-X", 
  "D09" = "D09-P", 
  "D11" = "D11-P", 
  "D12" = "D12-X", 
  "D13" = "D13-P", 
  "D15" = "D15-P", 
  "D16" = "D16-X"
)
combined[["Timepoint"]] <- Idents(combined)

# Save final merged Seurat object
saveRDS(combined, file = "combined.rds")

# DimPlots
dp.dataset <- DimPlot(
  object = combined, reduction = "umap", 
  group.by = "dataset", 
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Dataset of Origin"
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)
dp.timepoint <- DimPlot(
  object = combined, reduction = "umap", 
  group.by = "Timepoint", 
  cols = c("#CA4B0B", "#F8766D", "#DE8C00", "#D8C920", 
           "#A8C709", "#00BA38", "#18B9D5", "#2926D2", 
           "#C77CFF", "#FB61DE", "#7B5116", "#A5A2A9")
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Treatment Time"
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

png(
  filename = "../../Images/HumanReprogramming/xing+pana_UMAP_dataset.png", 
  width = 1420, height = 1000, res = 160
)
dp.dataset
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_UMAP_timepoint.png", 
  width = 1320, height = 1000, res = 160
)
dp.timepoint
dev.off()

# FeaturePlots
fp.nanog <- FeaturePlot(
  object = combined, 
  features = "NANOG", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.dppa4 <- FeaturePlot(
  object = combined, 
  features = "DPPA4", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.sema6a <- FeaturePlot(
  object = combined, 
  features = "SEMA6A", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.epcam <- FeaturePlot(
  object = combined, 
  features = "EPCAM", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

fp.col3a1 <- FeaturePlot(
  object = combined, 
  features = "COL3A1", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.mfap5 <- FeaturePlot(
  object = combined, 
  features = "MFAP5", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.tgfbi <- FeaturePlot(
  object = combined, 
  features = "TGFBI", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.mmp2 <- FeaturePlot(
  object = combined, 
  features = "MMP2", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

fp.cdkn1a <- FeaturePlot(
  object = combined, 
  features = "CDKN1A", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.cdkn2a <- FeaturePlot(
  object = combined, 
  features = "CDKN2A", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.tp53 <- FeaturePlot(
  object = combined, 
  features = "TP53", 
  order = TRUE
) + labs(
  x = "UMAP 1", y = "UMAP 2"
) + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_NANOG.png", 
  width = 1350, height = 1000, res = 140
)
fp.nanog
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_DPPA4.png", 
  width = 1300, height = 1000, res = 140
)
fp.dppa4
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_SEMA6A.png", 
  width = 1350, height = 1000, res = 140
)
fp.sema6a
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_EPCAM.png", 
  width = 1300, height = 1000, res = 140
)
fp.epcam
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_COL3A1.png", 
  width = 1300, height = 1000, res = 140
)
fp.col3a1
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_MFAP5.png", 
  width = 1300, height = 1000, res = 140
)
fp.mfap5
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_TGFBI.png", 
  width = 1300, height = 1000, res = 140
)
fp.tgfbi
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_MMP2.png", 
  width = 1300, height = 1000, res = 140
)
fp.mmp2
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_CDKN1A.png", 
  width = 1300, height = 1000, res = 140
)
fp.cdkn1a
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_CDKN2A.png", 
  width = 1300, height = 1000, res = 140
)
fp.cdkn2a
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_FeaturePlot_P53.png", 
  width = 1350, height = 1000, res = 140
)
fp.tp53
dev.off()


