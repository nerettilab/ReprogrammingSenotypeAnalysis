# Set working directory
setwd("Data/MergedHuman/")

# Load packages and Seurat object
library(dplyr)
library(stringr)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(monocle3)
library(shiny)

combined <- readRDS("combined.rds")
Idents(combined) <- "Timepoint"

# Create CellDataSet object to use with Monocle
cds.combined <- as.cell_data_set(combined, default.reduction = "umap")

cds.combined$ident <- NULL
fData(cds.combined) <- data.frame(
  gene_short_name = rownames(exprs(cds.combined)), 
  row.names = rownames(exprs(cds.combined))
)

cds.combined <- estimate_size_factors(cds = cds.combined)
cds.combined <- cluster_cells(cds.combined, reduction = "UMAP")

# Helper function to identify root principal point
get_earliest_principal_node <- function(
    cds = cds, 
    timepoint = "D00", 
    partition = "1"
) {
  cell_ids <- which(colData(cds)[, "Timepoint"] == timepoint & cds@clusters[["UMAP"]]$partitions == partition)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds),])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# Learn trajectory
cds.combined <- learn_graph(cds.combined, close_loop = FALSE)
cds.combined <- order_cells(cds.combined, root_pr_nodes = "Y_265")

saveRDS(cds.combined, file = "cds.combined.rds")

dp.pseudo <- plot_cells(
  cds = cds.combined, reduction_method = "UMAP", 
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.90, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, 
  graph_label_size = 3, cell_size = 0.35
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

png(
  filename = "../../Images/HumanReprogramming/xing+pana_pseudotime_full.png", 
  width = 1600, height = 1000, res = 200
)
dp.pseudo
dev.off()

## Compare gene expression between branches
# Somatic trajectory
cds.subset.som <- choose_graph_segments(
  cds = cds.combined, 
  starting_pr_node = "Y_34", 
  ending_pr_nodes = c("Y_322", "Y_292", "Y_268", "Y_75", 
                      "Y_2", "Y_20", "Y_31", "Y_405", 
                      "Y_166", "Y_369", "Y_366", "Y_62", 
                      "Y_73", "Y_337", "Y_473", "Y_119", 
                      "Y_118", "Y_448", "Y_188", "Y_221", 
                      "Y_211", "Y_143", "Y_170", "Y_401", 
                      "Y_193", "Y_194", "Y_527", "Y_471", 
                      "Y_524"), 
  clear_cds = FALSE
)

saveRDS(cds.subset.som, file = "cds.subset.som.rds")

cds.subset.som1 <- choose_graph_segments(
  cds = cds.combined, 
  starting_pr_node = "Y_34", 
  ending_pr_nodes = c("Y_322", "Y_292", "Y_268", 
                      "Y_75", "Y_2", "Y_20", 
                      "Y_31", "Y_405", "Y_524"), 
  clear_cds = FALSE
)

saveRDS(cds.subset.som1, file = "cds.subset.som1.rds")

cds.subset.som2 <- choose_graph_segments(
  cds = cds.combined, 
  starting_pr_node = "Y_34", 
  ending_pr_nodes = c("Y_322", "Y_292", "Y_268", "Y_75", 
                      "Y_2", "Y_20", "Y_31", "Y_405", 
                      "Y_166", "Y_369", "Y_366", "Y_62", 
                      "Y_73", "Y_337", "Y_473", "Y_119", 
                      "Y_118", "Y_448", "Y_188", "Y_221", 
                      "Y_211", "Y_143", "Y_170", "Y_401", 
                      "Y_193", "Y_194", "Y_527", "Y_471"), 
  clear_cds = FALSE
)

saveRDS(cds.subset.som2, file = "cds.subset.som2.rds")

# Pluripotent trajectory
cds.subset.plu <- choose_graph_segments(
  cds = cds.combined, 
  starting_pr_node = "Y_34", 
  ending_pr_nodes = c("Y_322", "Y_292", "Y_268", "Y_75", 
                      "Y_2", "Y_20", "Y_31", "Y_405", 
                      "Y_166", "Y_369", "Y_366", "Y_62", 
                      "Y_73", "Y_337", "Y_473", "Y_119", 
                      "Y_452", "Y_487", "Y_602", "Y_496", 
                      "Y_552", "Y_578", "Y_612"), 
  clear_cds = FALSE
)

saveRDS(cds.subset.plu, file = "cds.subset.plu.rds")

gc()

# Pseudotime plots of subsetted objects
dp.som.pseudotime <- plot_cells(
  cds = cds.subset.som, reduction_method = "UMAP", 
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.85, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, cell_size = 0.35
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

dp.som1.pseudotime <- plot_cells(
  cds = cds.subset.som1, reduction_method = "UMAP",
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.85, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, cell_size = 0.35
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

dp.som2.pseudotime <- plot_cells(
  cds = cds.subset.som2, reduction_method = "UMAP",
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.85, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, cell_size = 0.35
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

dp.plu.pseudotime <- plot_cells(
  cds = cds.subset.plu, reduction_method = "UMAP", 
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.85, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, cell_size = 0.35
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

png(
  filename = "../../Images/HumanReprogramming/xing+pana_pseudotime_somatic.png", 
  width = 900, height = 700, res = 120
)
dp.som.pseudotime
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_pseudotime_somatic1.png", 
  width = 780, height = 440, res = 120
)
dp.som1.pseudotime
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_pseudotime_somatic2.png", 
  width = 900, height = 700, res = 120
)
dp.som2.pseudotime
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_pseudotime_pluripotent.png", 
  width = 1070, height = 560, res = 120
)
dp.plu.pseudotime
dev.off()

# Plot genes in pseudotime
sen_genes <- c("CDKN1A", "CDKN2A", "TP53")
ylims <- list(
  CDKN1A = c(0.1, max(counts(cds.combined)["CDKN1A",])), 
  CDKN2A = c(0.1, max(counts(cds.combined)["CDKN2A",])), 
  TP53 = c(0.1, max(counts(cds.combined)["TP53",]))
)

kp.som <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.som[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#CA4B0B", "#F8766D", "#DE8C00", "#D8C920", 
               "#A8C709", "#00BA38", "#18B9D5", "#2926D2", 
               "#C77CFF", "#FB61DE", "#7B5116", "#A5A2A9")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})

kp.som1 <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.som1[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#CA4B0B", "#F8766D", "#DE8C00", "#D8C920", 
               "#A8C709", "#00BA38", "#18B9D5", "#2926D2", 
               "#C77CFF", "#FB61DE", "#7B5116", "#A5A2A9")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})

kp.som2 <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.som2[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#CA4B0B", "#F8766D", "#DE8C00", "#D8C920", 
               "#A8C709", "#00BA38", "#18B9D5", "#2926D2", 
               "#C77CFF", "#FB61DE", "#7B5116", "#A5A2A9")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})

kp.plu <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.plu[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#CA4B0B", "#F8766D", "#DE8C00", "#D8C920", 
               "#A8C709", "#00BA38", "#18B9D5", "#2926D2", 
               "#C77CFF", "#FB61DE", "#7B5116", "#A5A2A9")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})

png(
  filename = "../../Images/HumanReprogramming/xing+pana_kineticplot_somatic.png", 
  width = 996, height = 1000, res = 120
)
wrap_plots(
  kp.som, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_kineticplot_somatic1.png", 
  width = 770, height = 1000, res = 120
)
wrap_plots(
  kp.som1, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_kineticplot_somatic2.png", 
  width = 1100, height = 1000, res = 120
)
wrap_plots(
  kp.som2, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()
png(
  filename = "../../Images/HumanReprogramming/xing+pana_kineticplot_pluripotent.png", 
  width = 1100, height = 1000, res = 120
)
wrap_plots(
  kp.plu, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()


