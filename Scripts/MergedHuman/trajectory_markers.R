# Set working directory
setwd("Data/MergedHuman")

# Load packages and objects
library(Seurat)
library(presto)
library(ggplot2)
library(patchwork)
library(monocle3)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)

combined <- readRDS("combined.rds")
cds.subset.som1 <- readRDS("cds.subset.som1.rds")
cds.subset.som2 <- readRDS("cds.subset.som2.rds")
cds.subset.plu <- readRDS("cds.subset.plu.rds")


# Annotate which cells in Seurat object are part of which trajectory
cells.som1 <- colnames(cds.subset.som1)
cells.som2 <- colnames(cds.subset.som2)
cells.plu <- colnames(cds.subset.plu)
metadata <- combined@meta.data
metadata$Trajectory <- NA

metadata[cells.som1, "Trajectory"] <- "Somatic 1" 
metadata[cells.som2, "Trajectory"] <- "Somatic 2"
metadata[cells.plu, "Trajectory"] <- "Pluripotent"
metadata[intersect(cells.som2, cells.plu), "Trajectory"] <- "Before Split"
combined@meta.data <- metadata

dp.traj <- DimPlot(
  object = combined, reduction = "umap", 
  group.by = "Trajectory"
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Trajectory Branch"
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)
dp.traj.split <- DimPlot(
  object = combined, reduction = "umap", 
  group.by = "Trajectory", split.by = "dataset"
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Trajectory Branch"
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

png(
  filename = "../../Images/HumanReprogramming/xing+pana_UMAP_trajectorybranch.png", 
  width = 1450, height = 1000, res = 160
)
dp.traj
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_UMAP_trajectorybranch_splitbydataset.png", 
  width = 2390, height = 1000, res = 160
)
dp.traj.split
dev.off()

# Garbage collection
rm(cds.subset.plu, cds.subset.som1, cds.subset.som2, metadata, cells.som1, cells.som2, cells.plu)
gc(reset = TRUE)


## Violin plots of genes of interest
## Pluripotent: DPPA4, SEMA6A, EPCAM
## Somatic: MFAP5, TGFBI, MMP2
combined.aftersplit <- subset(combined, Trajectory != "Before Split")
timepoints.aftersplit <- c("D07-P", "D08-X", "D09-P", "D11-P", "D12-X", "D13-P", "D15-P", "D16-X")
genes.of.interest <- c("NANOG", "DPPA4", "SEMA6A", "EPCAM", "COL3A1", "MFAP5", "TGFBI", "MMP2", "CDKN1A", "CDKN2A", "TP53")

# Create data frames for statistical significance markers
sig.df.plu_vs_som1 <- data.frame(matrix(
  ncol = length(genes.of.interest), 
  nrow = length(timepoints.aftersplit)
))
colnames(sig.df.plu_vs_som1) <- genes.of.interest
rownames(sig.df.plu_vs_som1) <- timepoints.aftersplit
sig.df.plu_vs_som2 <- sig.df.plu_vs_som1
sig.df.som1_vs_som2 <- sig.df.plu_vs_som1


# Determine statistical significance of differential gene 
# expression between pluripotent and somatic trajectories at 
# each timepoint
Idents(combined.aftersplit) <- "Timepoint"
metadata.aftersplit <- combined.aftersplit@meta.data
for (timepoint in rownames(sig.df.plu_vs_som1)) {
  submeta <- metadata.aftersplit[which(metadata.aftersplit$Timepoint == timepoint),]
  submeta.plu <- nrow(submeta[which(submeta$Trajectory == "Pluripotent"),])
  submeta.som1 <- nrow(submeta[which(submeta$Trajectory == "Somatic 1"),])
  submeta.som2 <- nrow(submeta[which(submeta$Trajectory == "Somatic 2"),])
  
  # Perform Wilcoxon rank sum test if there are enough cells in each timepoint
  test.wilcox.plu_vs_som1 <- NA
  test.wilcox.plu_vs_som2 <- NA
  test.wilcox.som1_vs_som2 <- NA
  
  if (submeta.plu >= 3 & submeta.som1 >= 3) {
    test.wilcox.plu_vs_som1 <- FindMarkers(
      object = combined.aftersplit,
      group.by = "Trajectory", test.use = "wilcox",
      ident.1 = "Pluripotent", ident.2 = "Somatic 1",
      subset.ident = timepoint
    )
  }
  
  if (submeta.plu >= 3 & submeta.som2 >= 3) {
    test.wilcox.plu_vs_som2 <- FindMarkers(
      object = combined.aftersplit,
      group.by = "Trajectory", test.use = "wilcox",
      ident.1 = "Pluripotent", ident.2 = "Somatic 2",
      subset.ident = timepoint
    )
  }
  
  if (submeta.som1 >= 3 & submeta.som2 >= 3) {
    test.wilcox.som1_vs_som2 <- FindMarkers(
      object = combined.aftersplit,
      group.by = "Trajectory", test.use = "wilcox",
      ident.1 = "Somatic 1", ident.2 = "Somatic 2",
      subset.ident = timepoint
    )
  }
  
  
  # Assign asterisks based on p-value
  for (gene in colnames(sig.df.plu_vs_som1)) {
    if (gene %in% rownames(test.wilcox.plu_vs_som1)) {
      pval <- test.wilcox.plu_vs_som1[gene, "p_val_adj"]

      if (pval < 0.001) {
        sig.df.plu_vs_som1[timepoint, gene] <- "***"
      }
      else if (pval < 0.01) {
        sig.df.plu_vs_som1[timepoint, gene] <- "**"
      }
      else if (pval < 0.05) {
        sig.df.plu_vs_som1[timepoint, gene] <- "*"
      }
    }
  }
  
  for (gene in colnames(sig.df.plu_vs_som2)) {
    if (gene %in% rownames(test.wilcox.plu_vs_som2)) {
      pval <- test.wilcox.plu_vs_som2[gene, "p_val_adj"]
      
      if (pval < 0.001) {
        sig.df.plu_vs_som2[timepoint, gene] <- "***"
      }
      else if (pval < 0.01) {
        sig.df.plu_vs_som2[timepoint, gene] <- "**"
      }
      else if (pval < 0.05) {
        sig.df.plu_vs_som2[timepoint, gene] <- "*"
      }
    }
  }
  
  for (gene in colnames(sig.df.som1_vs_som2)) {
    if (gene %in% rownames(test.wilcox.som1_vs_som2)) {
      pval <- test.wilcox.som1_vs_som2[gene, "p_val_adj"]
      
      if (pval < 0.001) {
        sig.df.som1_vs_som2[timepoint, gene] <- "***"
      }
      else if (pval < 0.01) {
        sig.df.som1_vs_som2[timepoint, gene] <- "**"
      }
      else if (pval < 0.05) {
        sig.df.som1_vs_som2[timepoint, gene] <- "*"
      }
    }
  }
  

  gc()
}

sigtables <- list(sig.df.plu_vs_som1, sig.df.som1_vs_som2, sig.df.plu_vs_som2)


# Determine differential expression of genes at end of
# trajectory compared to at point of trajectory 
# divergence for each trajectory branch
Idents(combined.aftersplit) <- "Trajectory"
test.wilcox.plu <- FindMarkers(
  object = combined.aftersplit,
  group.by = "Timepoint",
  ident.1 = "D16-X", ident.2 = "D07-P",
  subset.ident = "Pluripotent"
)
test.wilcox.som1 <- FindMarkers(
  object = combined.aftersplit,
  group.by = "Timepoint",
  ident.1 = "D16-X", ident.2 = "D07-P",
  subset.ident = "Somatic 1"
)
test.wilcox.som2 <- FindMarkers(
  object = combined.aftersplit,
  group.by = "Timepoint",
  ident.1 = "D16-X", ident.2 = "D07-P",
  subset.ident = "Somatic 2"
)

gc()

# Create violin plots for time points after trajectory split
Idents(combined.aftersplit) <- "Timepoint"
for (gene in colnames(sig.df.plu_vs_som1)) {
  y.max <- max(FetchData(combined.aftersplit, vars = gene)) # Determine bounds of y-axis
  
  
  # Create initial violin plots
  vp <- VlnPlot(
    object = combined.aftersplit, features = gene, 
    cols = c("#7CAE00", "#00BFC4", "#C77CFF"), 
    pt.size = 0.001, idents = timepoints.aftersplit, 
    group.by = "Timepoint", split.by = "Trajectory"
  ) + scale_x_discrete(drop = TRUE) + ylim(0, y.max * 1.34) + theme(
    text = element_text(size = 30), 
    axis.text = element_text(size = 21), 
    legend.key.size = unit(1.15, "cm")
  ) + labs(x = "Timepoint")
  
  
  # Variables to handle spacing of statistical significance markers
  sigbar.counter <- y.max * 1.06
  sigbar.increment <- sigbar.counter * 0.022
  sigbar.max.y <- sigbar.counter
  
  # Add statistical significance markers between 
  # trajectory branches at each timepoint
  for (i in seq_along(timepoints.aftersplit)) {
    sigbar.counter <- y.max * 1.06
    
    for (sigtable in seq_along(sigtables)) {
      if (!is.na(sigtables[[sigtable]][i, gene])) {
        # Determine horizontal boundaries of statistical significance markers
        xstart = switch(
          sigtable, 
          "1" = (i - 0.3), 
          "2" = i, 
          "3" = (i - 0.3)
        )
        xend = switch(
          sigtable,
          "1" = xstart + 0.3,
          "2" = xstart + 0.3,
          "3" = xstart + 0.6
        )
        
        # Add statistical significance markers to violin plot
        vp <- vp + annotate(
          "segment", 
          x = xstart , xend = xend, 
          y = sigbar.counter, yend = sigbar.counter
        ) + annotate(
          "text", 
          x = (xstart + xend) / 2, 
          y = sigbar.counter + 0.75*sigbar.increment, 
          label = sigtables[[sigtable]][i, gene]
        )
        
        sigbar.counter <- sigbar.counter + 2*sigbar.increment
        if (sigbar.counter > sigbar.max.y) {
          sigbar.max.y <- sigbar.counter
        }
      }
    }
  }
  
  
  # Add statistical significance markers over time for each 
  # trajectory branch
  sig.time.plu <- NA
  sig.time.som1 <- NA
  sig.time.som2 <- NA
  
  if (gene %in% rownames(test.wilcox.plu)) {
    pval <- test.wilcox.plu[gene, "p_val_adj"]
    if (pval < 0.001) { sig.time.plu <- "***" }
    else if (pval < 0.01) { sig.time.plu <- "**" }
    else if (pval < 0.05) { sig.time.plu <- "*" }
  }
  if (gene %in% rownames(test.wilcox.som1)) {
    pval <- test.wilcox.som1[gene, "p_val_adj"]
    if (pval < 0.001) { sig.time.som1 <- "***" }
    else if (pval < 0.01) { sig.time.som1 <- "**" }
    else if (pval < 0.05) { sig.time.som1 <- "*" }
  }
  if (gene %in% rownames(test.wilcox.som2)) {
    pval <- test.wilcox.som2[gene, "p_val_adj"]
    if (pval < 0.001) { sig.time.som2 <- "***" }
    else if (pval < 0.01) { sig.time.som2 <- "**" }
    else if (pval < 0.05) { sig.time.som2 <- "*" }
  }
  
  sigtimes <- c(sig.time.plu, sig.time.som1, sig.time.som2)
  sigbar.counter <- sigbar.max.y
  
  for (sigtime in seq_along(sigtimes)) {
    if (!is.na(sigtimes[sigtime])) {
      seq.begin <- seq_along(timepoints.aftersplit)[1] - 0.3
      seq.end <- seq_along(timepoints.aftersplit)[length(seq_along(timepoints.aftersplit))] - 0.3
      
      vp <- vp + annotate(
        "segment",
        x = seq.begin + 0.3*(sigtime - 1), 
        xend = seq.end + 0.3*(sigtime - 1), 
        y = sigbar.counter, yend = sigbar.counter
      ) + annotate(
        "text", 
        x = ((seq.begin + 0.3*(sigtime - 1)) + (seq.end + 0.3*(sigtime - 1))) / 2, 
        y = sigbar.counter + 0.75*sigbar.increment, 
        label = sigtimes[sigtime]
      )
      
      sigbar.counter <- sigbar.counter + 2*sigbar.increment
    }
  }
  
  
  # Assign variable to completed violin plot
  assign(paste0("vp.", gene), vp)
  
  gc()
}

# Save images of violin plots
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_NANOG.png", 
  width = 1550, height = 1100, res = 140
)
vp.NANOG
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_DPPA4.png", 
  width = 1550, height = 1100, res = 140
)
vp.DPPA4
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_SEMA6A.png", 
  width = 1550, height = 1100, res = 140
)
vp.SEMA6A
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_EPCAM.png", 
  width = 1550, height = 1100, res = 140
)
vp.EPCAM
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_COL3A1.png", 
  width = 1550, height = 1100, res = 140
)
vp.COL3A1
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_MFAP5.png", 
  width = 1550, height = 1100, res = 140
)
vp.MFAP5
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_TGFBI.png", 
  width = 1550, height = 1100, res = 140
)
vp.TGFBI
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_MMP2.png", 
  width = 1550, height = 1100, res = 140
)
vp.MMP2
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_CDKN1A.png", 
  width = 1550, height = 1100, res = 140
)
vp.CDKN1A
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_CDKN2A.png", 
  width = 1550, height = 1100, res = 140
)
vp.CDKN2A
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_ViolinPlot_P53.png", 
  width = 1550, height = 1100, res = 140
)
vp.TP53
dev.off()


# Identify genes upregulated in pluripotent/somatic branches compared to start of time-course (Day 0)
markers.plu <- FindMarkers(
  object = combined[["RNA"]],
  cells.1 = Cells(subset(combined, Trajectory == "Pluripotent")), 
  cells.2 = Cells(subset(combined, Timepoint == "D00" & Trajectory == "Before Split")), 
  logfc.threshold = 0.5, min.pct = 0.30, only.pos = TRUE
)
markers.plu <- subset(markers.plu, p_val_adj < 0.05)
gc()

markers.som <- FindMarkers(
  object = combined[["RNA"]],
  cells.1 = Cells(subset(combined, Trajectory == "Somatic 1" | Trajectory == "Somatic 2")), 
  cells.2 = Cells(subset(combined, Timepoint == "D00" & Trajectory == "Before Split")), 
  logfc.threshold = 0.5, min.pct = 0.30, only.pos = TRUE
)
markers.som <- subset(markers.som, p_val_adj < 0.05)
gc()

markers.som1 <- FindMarkers(
  object = combined[["RNA"]],
  cells.1 = Cells(subset(combined, Trajectory == "Somatic 1")), 
  cells.2 = Cells(subset(combined, Timepoint == "D00" & Trajectory == "Before Split")), 
  logfc.threshold = 0.5, min.pct = 0.30, only.pos = TRUE
)
markers.som1 <- subset(markers.som1, p_val_adj < 0.05)
gc()

markers.som2 <- FindMarkers(
  object = combined[["RNA"]],
  cells.1 = Cells(subset(combined, Trajectory == "Somatic 2")), 
  cells.2 = Cells(subset(combined, Timepoint == "D00" & Trajectory == "Before Split")), 
  logfc.threshold = 0.5, min.pct = 0.30, only.pos = TRUE
)
markers.som2 <- subset(markers.som2, p_val_adj < 0.05)
gc()

# Convert marker gene names to Entrez IDs
genedf.plu <- bitr(
  rownames(markers.plu), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)
genedf.som <- bitr(
  rownames(markers.som), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)
genedf.som1 <- bitr(
  rownames(markers.som1), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)
genedf.som2 <- bitr(
  rownames(markers.som2), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

# Obtain background gene set
bg.genes <- bitr(
  Features(combined, layer = "data"), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

# GO analysis
go.plu <- enrichGO(
  gene = as.character(genedf.plu$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)
go.som <- enrichGO(
  gene = as.character(genedf.som$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)
go.som1 <- enrichGO(
  gene = as.character(genedf.som1$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)
go.som2 <- enrichGO(
  gene = as.character(genedf.som2$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)

# Plots of GO analysis
dp.go.plu <- dotplot(go.plu, showCategory = 20) + ggtitle("GO Biological Processes (Pluripotent)")
dp.go.som <- dotplot(go.som, showCategory = 20) + ggtitle("GO Biological Processes (Somatic)")
dp.go.som1 <- dotplot(go.som1, showCategory = 20) + ggtitle("GO Biological Processes (Somatic 1)")
dp.go.som2 <- dotplot(go.som2, showCategory = 20) + ggtitle("GO Biological Processes (Somatic 2)")

png(
  filename = "../../Images/HumanReprogramming/human_trajectory_GO.png", 
  width = 2050, height = 2050, res = 160
)
dp.go.plu | dp.go.som
dev.off()

png(
  filename = "../../Images/HumanReprogramming/human_trajectory_GO_somaticclusters.png", 
  width = 2050, height = 2050, res = 160
)
dp.go.som1 | dp.go.som2
dev.off()

# Identify which enriched GO processes p53 is a part of
res.go.plu <- go.plu@result
for (i in 1:20) {
  res.list <- str_split(res.go.plu[i,"geneID"], pattern = "/")[[1]]
  if ("TP53" %in% res.list) {
    print(res.go.plu[i, "Description"])
  }
}


