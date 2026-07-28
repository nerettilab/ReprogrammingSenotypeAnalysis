# Set working directory
setwd("Data/Schiebinger2019/")

# Load packages and objects
library(Seurat)
library(presto)
library(ggplot2)
library(patchwork)
library(monocle3)
library(org.Mm.eg.db)
library(clusterProfiler)
library(stringr)

scbg <- readRDS("scbg_final.rds")
cds.subset.stro <- readRDS("cds.subset.stro.rds")
cds.subset.neur <- readRDS("cds.subset.neur.rds")
cds.subset.trop <- readRDS("cds.subset.trop.rds")
cds.subset.plu <- readRDS("cds.subset.plu.rds")

# Annotate which cells in Seurat object are part of which trajectory/lineage
cells.stro <- colnames(cds.subset.stro)
cells.neur <- colnames(cds.subset.neur)
cells.trop <- colnames(cds.subset.trop)
cells.plu <- colnames(cds.subset.plu)
metadata <- scbg@meta.data
metadata$Trajectory <- NA

metadata[cells.stro, "Trajectory"] <- "Stromal" 
metadata[cells.neur, "Trajectory"] <- "Neural"
metadata[cells.trop, "Trajectory"] <- "Trophoblast"
metadata[cells.plu, "Trajectory"] <- "Pluripotent"
metadata[intersect(cells.trop, cells.plu), "Trajectory"] <- "Before Split"
metadata[setdiff(rownames(metadata[which(metadata$Trajectory == "Before Split"),]), cells.stro), "Trajectory"] <- "MET"
scbg@meta.data <- metadata

scbg$Trajectory <- factor(
  scbg$Trajectory, 
  levels = c("Before Split", "Stromal", "MET", "Neural", "Trophoblast", "Pluripotent")
)

dp.traj <- DimPlot(
  object = scbg, reduction = "umap", 
  cols = c("#F8766D", "#B79F00", "#00BA38", "#18B9D5", "#2926D2", "#C77CFF"), 
  group.by = "Trajectory", raster = TRUE
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
  filename = "../../Images/MouseReprogramming/scbg_UMAP_trajectorybranch.png", 
  width = 1080, height = 750, res = 120
)
dp.traj
dev.off()

rm(cds.subset.stro, cds.subset.neur, cds.subset.trop, cds.subset.plu, metadata, cells.stro, cells.neur, cells.trop, cells.plu)
gc(reset = TRUE)


## Violin plots of genes of interest
## Pluripotent: Dppa5a, Tdgf1, Dusp9,
## Stromal: Mmp2, Fmod, Cd47

# Group MET and Pluripotent cells together
scbg.plu_vs_stro <- subset(scbg, Trajectory == "Stromal" | Trajectory == "MET" | Trajectory == "Pluripotent")
scbg.plu_vs_stro$Trajectory <- as.character(scbg.plu_vs_stro$Trajectory)

metadata.plu_vs_stro <- scbg.plu_vs_stro@meta.data
metadata.plu_vs_stro[which(metadata.plu_vs_stro$Trajectory == "MET" | metadata.plu_vs_stro$Trajectory == "Pluripotent"), "Trajectory"] <- "MET/Pluripotent"
scbg.plu_vs_stro@meta.data <- metadata.plu_vs_stro

# Create data frame for statistical significance markers
timepoints.plu_vs_stro <- c(
  "D08", "D09", "D10", "D11", "D12", "D13", 
  "D14", "D15", "D16", "D17", "D18"
)
genes.of.interest <- c("Nanog", "Dppa5a", "Tdgf1", "Dusp9", "Col3a1", "Mmp2", "Fmod", "Cd47", "Cdkn1a", "Cdkn2a", "Trp53")

sig.df <- data.frame(matrix(
  ncol = length(genes.of.interest), 
  nrow = length(timepoints.plu_vs_stro)
))
colnames(sig.df) <- genes.of.interest
rownames(sig.df) <- timepoints.plu_vs_stro

# Determine statistical significance of differential gene 
# expression between trajectory branches at each 
# timepoint
Idents(scbg.plu_vs_stro) <- "Timepoint"
for (timepoint in rownames(sig.df)) {
  submeta <- metadata.plu_vs_stro[which(metadata.plu_vs_stro$Timepoint == timepoint),]
  submeta.plu <- nrow(submeta[which(submeta$Trajectory == "MET/Pluripotent"),])
  submeta.stro <- nrow(submeta[which(submeta$Trajectory == "Stromal"),])
  
  # Perform Wilcoxon rank sum test if there are enough cells in each timepoint
  test.wilcox <- NA
  if (submeta.plu >= 3 & submeta.stro >= 3) {
    test.wilcox <- FindMarkers(
      object = scbg.plu_vs_stro,
      group.by = "Trajectory", test.use = "wilcox",
      ident.1 = "MET/Pluripotent", ident.2 = "Stromal",
      subset.ident = timepoint
    )
  }
  
  # Assign asterisks based on p-value
  for (gene in colnames(sig.df)) {
    if (gene %in% rownames(test.wilcox)) {
      pval <- test.wilcox[gene, "p_val_adj"]
      
      if (pval < 0.001) {
        sig.df[timepoint, gene] <- "***"
      }
      else if (pval < 0.01) {
        sig.df[timepoint, gene] <- "**"
      }
      else if (pval < 0.05) {
        sig.df[timepoint, gene] <- "*"
      }
    }
  }
  
  gc()
}


# Determine differential expression of genes at end of
# trajectory compared to at point of trajectory 
# divergence for each trajectory branch
Idents(scbg.plu_vs_stro) <- "Trajectory"
test.wilcox.plu <- FindMarkers(
  object = scbg.plu_vs_stro,
  group.by = "Timepoint",
  ident.1 = "D18", ident.2 = "D08",
  subset.ident = "MET/Pluripotent"
)
test.wilcox.stro <- FindMarkers(
  object = scbg.plu_vs_stro,
  group.by = "Timepoint",
  ident.1 = "D18", ident.2 = "D08",
  subset.ident = "Stromal"
)

gc()

# Create violin plots
Idents(scbg.plu_vs_stro) <- "Timepoint"
for (gene in colnames(sig.df)) {
  y.max <- max(FetchData(scbg.plu_vs_stro, vars = gene)) # Determine bounds of y-axis
  
  # Create initial violin plot
  vp <- VlnPlot(
    object = scbg.plu_vs_stro, features = gene, 
    cols = c("#C77CFF", "#B79F00"), 
    pt.size = 0.001, idents = timepoints.plu_vs_stro, 
    group.by = "Timepoint", split.by = "Trajectory", 
    raster = TRUE
  ) + scale_x_discrete(drop = TRUE) + ylim(0, y.max * 1.21) + theme(
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
  for (i in seq_along(timepoints.plu_vs_stro)) {
    if (!is.na(sig.df[i, gene])) {
      sigbar.counter <- y.max * 1.06
      
      vp <- vp + annotate(
        "segment", 
        x = i - 0.21 , xend = i + 0.21, 
        y = sigbar.counter, yend = sigbar.counter
      ) + annotate(
        "text", 
        x = i, 
        y = sigbar.counter + 0.75*sigbar.increment, 
        label = sig.df[i, gene]
      )
      
      sigbar.counter <- sigbar.counter + 2*sigbar.increment
      if (sigbar.counter > sigbar.max.y) {
        sigbar.max.y <- sigbar.counter
      }
    }
  }
  
  
  # Add statistical significance markers over time for each 
  # trajectory branch
  sig.time.plu <- NA
  sig.time.stro <- NA
  
  if (gene %in% rownames(test.wilcox.plu)) {
    pval <- test.wilcox.plu[gene, "p_val_adj"]
    if (pval < 0.001) { sig.time.plu <- "***" }
    else if (pval < 0.01) { sig.time.plu <- "**" }
    else if (pval < 0.05) { sig.time.plu <- "*" }
  }
  if (gene %in% rownames(test.wilcox.stro)) {
    pval <- test.wilcox.stro[gene, "p_val_adj"]
    if (pval < 0.001) { sig.time.stro <- "***" }
    else if (pval < 0.01) { sig.time.stro <- "**" }
    else if (pval < 0.05) { sig.time.stro <- "*" }
  }
  
  sigtimes <- c(sig.time.plu, sig.time.stro)
  sigbar.counter <- sigbar.max.y
  
  for (sigtime in seq_along(sigtimes)) {
    if (!is.na(sigtimes[sigtime])) {
      seq.begin <- seq_along(timepoints.plu_vs_stro)[1] - 0.21
      seq.end <- seq_along(timepoints.plu_vs_stro)[length(seq_along(timepoints.plu_vs_stro))] - 0.21
      
      vp <- vp + annotate(
        "segment",
        x = seq.begin + 0.42*(sigtime - 1), 
        xend = seq.end + 0.42*(sigtime - 1), 
        y = sigbar.counter, yend = sigbar.counter
      ) + annotate(
        "text", 
        x = ((seq.begin + 0.42*(sigtime - 1)) + (seq.end + 0.42*(sigtime - 1))) / 2, 
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
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_NANOG.png", 
  width = 1900, height = 1220, res = 140
)
vp.Nanog
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_DPPA5A.png", 
  width = 1900, height = 1220, res = 140
)
vp.Dppa5a
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_TDGF1.png", 
  width = 1900, height = 1220, res = 140
)
vp.Tdgf1
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_DUSP9.png", 
  width = 1900, height = 1220, res = 140
)
vp.Dusp9
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_COL3A1.png", 
  width = 1900, height = 1220, res = 140
)
vp.Col3a1
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_MMP2.png", 
  width = 1900, height = 1220, res = 140
)
vp.Mmp2
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_FMOD.png", 
  width = 1900, height = 1220, res = 140
)
vp.Fmod
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_CD47.png", 
  width = 1900, height = 1220, res = 140
)
vp.Cd47
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_CDKN1A.png", 
  width = 1900, height = 1220, res = 140
)
vp.Cdkn1a
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_CDKN2A.png", 
  width = 1900, height = 1220, res = 140
)
vp.Cdkn2a
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_ViolinPlot_P53.png", 
  width = 1900, height = 1220, res = 140
)
vp.Trp53
dev.off()

# Garbage collection
rm(sig.df, test.wilcox, test.wilcox.plu, test.wilcox.stro, vp, gene, genes.of.interest, i, pval, seq.begin, seq.end, sigtime, sigtimes, sig.time.plu, sig.time.stro, timepoint, y.max, sigbar.counter, sigbar.increment, sigbar.max.y, submeta, submeta.plu, submeta.stro, timepoints.plu_vs_stro, scbg.plu_vs_stro, metadata.plu_vs_stro)
gc(reset = TRUE)


# Identify markers of somatic and pluripotent trajectory branches
markers.plu <- FindMarkers(
  object = scbg[["RNA"]],
  cells.1 = Cells(subset(scbg, Trajectory == "Pluripotent")), 
  cells.2 = Cells(subset(scbg, Timepoint == "D00" & Trajectory == "Before Split")), 
  logfc.threshold = 0.5, min.pct = 0.30, only.pos = TRUE
)
markers.plu <- subset(markers.plu, p_val_adj < 0.05)
gc()

markers.stro <- FindMarkers(
  object = scbg[["RNA"]],
  cells.1 = Cells(subset(scbg, Trajectory == "Stromal")), 
  cells.2 = Cells(subset(scbg, Timepoint == "D00" & Trajectory == "Before Split")), 
  logfc.threshold = 0.5, min.pct = 0.30, only.pos = TRUE
)
markers.stro <- subset(markers.stro, p_val_adj < 0.05)
gc()

markers.som <- FindMarkers(
  object = scbg[["RNA"]],
  cells.1 = Cells(subset(scbg, Trajectory == "Stromal" | Trajectory == "Neural")), 
  cells.2 = Cells(subset(scbg, Timepoint == "D00" & Trajectory == "Before Split")), 
  logfc.threshold = 0.5, min.pct = 0.30, only.pos = TRUE
)
markers.som <- subset(markers.som, p_val_adj < 0.05)
gc()

# Convert marker gene names to Entrez IDs
genedf.plu <- bitr(
  rownames(markers.plu), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Mm.eg.db
)
genedf.stro <- bitr(
  rownames(markers.stro), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Mm.eg.db
)
genedf.som <- bitr(
  rownames(markers.som), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Mm.eg.db
)

# Obtain background gene set
bg.genes <- bitr(
  Features(scbg, layer = "data"), 
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = org.Mm.eg.db
)

# GO/KEGG analysis
go.plu <- enrichGO(
  gene = as.character(genedf.plu$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Mm.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)
go.stro <- enrichGO(
  gene = as.character(genedf.stro$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Mm.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)
go.som <- enrichGO(
  gene = as.character(genedf.som$ENTREZID), 
  universe = as.character(bg.genes$ENTREZID), 
  OrgDb = org.Mm.eg.db, 
  ont = "BP", # Biological process
  qvalueCutoff = 0.05, 
  readable = TRUE
)

# Plots of GO analysis
dp.go.plu <- dotplot(go.plu, showCategory = 20) + ggtitle("GO Biological Processes (Pluripotent)")
dp.go.stro <- dotplot(go.stro, showCategory = 20) + ggtitle("GO Biological Processes (Stromal)")
dp.go.som <- dotplot(go.som, showCategory = 20) + ggtitle("GO Biological Processes (Somatic)")

png(
  filename = "../../Images/MouseReprogramming/mouse_trajectory_GO_pluripotent_vs_stromal.png", 
  width = 2050, height = 1800, res = 160
)
dp.go.plu | dp.go.stro
dev.off()

png(
  filename = "../../Images/MouseReprogramming/mouse_trajectory_GO.png", 
  width = 2050, height = 1800, res = 160
)
dp.go.plu | dp.go.som
dev.off()

# Identify which enriched GO processes p53 is a part of
res.go.plu <- go.plu@result
for (i in 1:20) {
  res.list <- str_split(res.go.plu[i,"geneID"], pattern = "/")[[1]]
  if ("Trp53" %in% res.list) {
    print(res.go.plu[i, "Description"])
  }
}


