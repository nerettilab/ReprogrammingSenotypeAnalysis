# Set working directory
setwd("Data/Schiebinger2019/")

# Load packages and Seurat object
library(dplyr)
library(stringr)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(monocle3)
library(shiny)

scbg <- readRDS("scbg_final.rds")
Idents(scbg) <- "Timepoint"

# Create CellDataSet object to use with Monocle
cds.scbg <- as.cell_data_set(scbg, default.reduction = "umap")

cds.scbg$ident <- NULL
fData(cds.scbg) <- data.frame(
  gene_short_name = rownames(exprs(cds.scbg)), 
  row.names = rownames(exprs(cds.scbg))
)

cds.scbg <- estimate_size_factors(cds = cds.scbg)
cds.scbg <- cluster_cells(cds.scbg, reduction = "UMAP")

rm(scbg)
gc(reset = TRUE)

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
cds.scbg <- learn_graph(cds.scbg)
cds.scbg <- order_cells(cds.scbg, root_pr_nodes = "Y_4")

saveRDS(cds.scbg, file = "cds.scbg.rds")

dp.pseudo <- plot_cells(
  cds = cds.scbg, reduction_method = "UMAP", 
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.90, 
  label_branch_points = FALSE, 
  label_leaves = FALSE, 
  label_cell_groups = FALSE, 
  graph_label_size = 3, cell_size = 0.25
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

png(
  filename = "../../Images/MouseReprogramming/scbg_pseudotime_full.png", 
  width = 1600, height = 1000, res = 180
)
dp.pseudo
dev.off()


## Compare gene expression between branches
# Somatic trajectory
cds.subset.som <- choose_graph_segments(
  cds = cds.scbg, 
  starting_pr_node = "Y_4", 
  ending_pr_nodes = c("Y_100", "Y_1", "Y_72", "Y_106", 
                      "Y_905", "Y_915", "Y_885", "Y_935", 
                      "Y_1002", "Y_1064", "Y_927", "Y_1142", 
                      "Y_1059", "Y_1049", "Y_1063", "Y_1078", 
                      "Y_991", "Y_1219", "Y_1181", "Y_1162", 
                      "Y_1334", "Y_1170", "Y_1389", "Y_1339", 
                      "Y_1150", "Y_1252", "Y_1255", "Y_1366", 
                      "Y_1362", "Y_1217", "Y_1693", "Y_986", 
                      "Y_1240", "Y_1355", "Y_1343", "Y_1391", 
                      "Y_1484", "Y_1516", "Y_1508", "Y_1475", 
                      "Y_1533", "Y_1632", "Y_1569", "Y_1634", 
                      "Y_1631",  "Y_1804", "Y_332", "Y_848", 
                      "Y_873", "Y_1437", "Y_1768", "Y_1652", 
                      "Y_1670", "Y_1681", "Y_1678", "Y_1803", 
                      "Y_1774", "Y_1786", "Y_126", "Y_247", 
                      "Y_256", "Y_264", "Y_396", "Y_733", 
                      "Y_661", "Y_317", "Y_596", "Y_139", 
                      "Y_175", "Y_233", "Y_576", "Y_255", 
                      "Y_492", "Y_155", "Y_455", "Y_368", 
                      "Y_837"), 
  clear_cds = FALSE
)
saveRDS(cds.subset.som, file = "cds.subset.som.rds")
gc()

# Stromal lineage
cds.subset.stro <- choose_graph_segments(
  cds = cds.scbg, 
  starting_pr_node = "Y_4", 
  ending_pr_nodes = c("Y_100", "Y_1", "Y_72", "Y_106", 
                      "Y_905", "Y_915", "Y_885", "Y_935", 
                      "Y_1002", "Y_1064", "Y_927", "Y_1142", 
                      "Y_1059", "Y_1049", "Y_1063", "Y_1078", 
                      "Y_991", "Y_1219", "Y_1181", "Y_1162", 
                      "Y_1334", "Y_1170", "Y_1389", "Y_1339", 
                      "Y_1150", "Y_1252", "Y_1255", "Y_1366", 
                      "Y_1362", "Y_1217", "Y_1693", "Y_986", 
                      "Y_1240", "Y_1355", "Y_1343", "Y_1391", 
                      "Y_1484", "Y_1516", "Y_1508", "Y_1475", 
                      "Y_1533", "Y_1632", "Y_1569", "Y_1634", 
                      "Y_1631",  "Y_1437", "Y_1768", "Y_1652", 
                      "Y_1670", "Y_1681", "Y_1678", "Y_1803", 
                      "Y_1774", "Y_1786", "Y_126", "Y_247", 
                      "Y_256", "Y_264", "Y_396", "Y_733", 
                      "Y_661", "Y_317", "Y_596", "Y_139", 
                      "Y_175", "Y_233", "Y_576", "Y_255", 
                      "Y_492", "Y_155", "Y_455", "Y_368", 
                      "Y_837"), 
  clear_cds = FALSE
)
saveRDS(cds.subset.stro, file = "cds.subset.stro.rds")
gc()

# Neural lineage
cds.subset.neur <- choose_graph_segments(
  cds = cds.scbg, 
  starting_pr_node = "Y_4", 
  ending_pr_nodes = c("Y_100", "Y_1", "Y_72", "Y_106", 
                      "Y_905", "Y_915", "Y_885", "Y_935", 
                      "Y_1002", "Y_1064", "Y_927", "Y_1142", 
                      "Y_1059", "Y_1049", "Y_1063", "Y_1078", 
                      "Y_991", "Y_1219", "Y_1181", "Y_1162", 
                      "Y_1334", "Y_1170", "Y_1389", "Y_1339", 
                      "Y_1150", "Y_1252", "Y_1255", "Y_1366", 
                      "Y_1362", "Y_1217", "Y_1693", "Y_986", 
                      "Y_1240", "Y_1355", "Y_1343", "Y_1391", 
                      "Y_1484", "Y_1516", "Y_1508", "Y_1475", 
                      "Y_1533", "Y_1632", "Y_1569", "Y_1634", 
                      "Y_1631", "Y_1804", "Y_332", "Y_848", 
                      "Y_873"), 
  clear_cds = FALSE
)
saveRDS(cds.subset.neur, file = "cds.subset.neur.rds")
gc()

# Trophoblast lineage
cds.subset.trop <- choose_graph_segments(
  cds = cds.scbg, 
  starting_pr_node = "Y_4", 
  ending_pr_nodes = c("Y_100", "Y_1", "Y_72", "Y_106", 
                      "Y_905", "Y_915", "Y_885", "Y_935", 
                      "Y_1002", "Y_1064", "Y_927", "Y_1142", 
                      "Y_1059", "Y_1049", "Y_1063", "Y_1078", 
                      "Y_991", "Y_1219", "Y_1181", "Y_1162", 
                      "Y_1334", "Y_1170", "Y_1389", "Y_1339", 
                      "Y_1150", "Y_1252", "Y_1255", "Y_1366", 
                      "Y_1362", "Y_1217", "Y_1693", "Y_986", 
                      "Y_1240", "Y_1355", "Y_1343", "Y_1391", 
                      "Y_1484", "Y_1516", "Y_1508", "Y_1475", 
                      "Y_1533", "Y_1632", "Y_1569", "Y_1634", 
                      "Y_1631", "Y_1804", "Y_332",  "Y_674", 
                      "Y_862", "Y_287", "Y_153"), 
  clear_cds = FALSE
)
saveRDS(cds.subset.trop, file = "cds.subset.trop.rds")
gc()

# Pluripotent trajectory
cds.subset.plu <- choose_graph_segments(
  cds = cds.scbg, 
  starting_pr_node = "Y_4", 
  ending_pr_nodes = c("Y_100", "Y_1", "Y_72", "Y_106", 
                      "Y_905", "Y_915", "Y_885", "Y_935", 
                      "Y_1002", "Y_1064", "Y_927", "Y_1142", 
                      "Y_1059", "Y_1049", "Y_1063", "Y_1078", 
                      "Y_991", "Y_1219", "Y_1181", "Y_1162", 
                      "Y_1334", "Y_1170", "Y_1389", "Y_1339", 
                      "Y_1150", "Y_1252", "Y_1255", "Y_1366", 
                      "Y_1362", "Y_1217", "Y_1693", "Y_986", 
                      "Y_1240", "Y_1355", "Y_1343", "Y_1391", 
                      "Y_1484", "Y_1516", "Y_1508", "Y_1475", 
                      "Y_1533", "Y_1632", "Y_1569", "Y_1634", 
                      "Y_1631", "Y_1804", "Y_332", "Y_674", 
                      "Y_1939", "Y_791", "Y_656", "Y_403", 
                      "Y_374", "Y_1830", "Y_273", "Y_1860", 
                      "Y_1868", "Y_378", "Y_754", "Y_654", 
                      "Y_656", "Y_402", "Y_637", "Y_499", 
                      "Y_504", "Y_612", "Y_727", "Y_503", 
                      "Y_763", "Y_726", "Y_548", "Y_572", 
                      "Y_419"), 
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
  label_cell_groups = FALSE, cell_size = 0.25
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)
dp.stro.pseudotime <- plot_cells(
  cds = cds.subset.stro, reduction_method = "UMAP", 
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.85, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, cell_size = 0.25
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)
dp.neur.pseudotime <- plot_cells(
  cds = cds.subset.neur, reduction_method = "UMAP", 
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.85, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, cell_size = 0.25
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)
dp.trop.pseudotime <- plot_cells(
  cds = cds.subset.trop, reduction_method = "UMAP", 
  color_cells_by = "pseudotime", 
  show_trajectory_graph = TRUE, 
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 0.85, 
  label_branch_points = FALSE, label_leaves = FALSE, 
  label_cell_groups = FALSE, cell_size = 0.25
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
  label_cell_groups = FALSE, cell_size = 0.25
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

png(
  filename = "../../Images/MouseReprogramming/scbg_pseudotime_som.png", 
  width = 880, height = 750, res = 120
)
dp.som.pseudotime
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_pseudotime_stro.png", 
  width = 850, height = 700, res = 120
)
dp.stro.pseudotime
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_pseudotime_neur.png", 
  width = 880, height = 750, res = 120
)
dp.neur.pseudotime
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_pseudotime_trop.png", 
  width = 1020, height = 630, res = 120
)
dp.trop.pseudotime
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_pseudotime_plu.png", 
  width = 1170, height = 580, res = 120
)
dp.plu.pseudotime
dev.off()

gc()


## Plot genes in pseudotime
sen_genes <- c("Cdkn1a", "Cdkn2a", "Trp53")
ylims <- list(
  Cdkn1a = c(0.1, max(counts(cds.scbg)["Cdkn1a",])), 
  Cdkn2a = c(0.1, max(counts(cds.scbg)["Cdkn2a",])), 
  Trp53 = c(0.1, max(counts(cds.scbg)["Trp53",]))
)

kp.som <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.som[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#8B0000", "#CA4B2C", "#F8766D", "#DE8C00", 
               "#D8C920", "#A8C709", "#094500", "#00BA38", 
               "#1DDED0", "#224AD9", "#000168", "#6E15D4", 
               "#C77CFF", "#FB61DE", "#C0BC87", "#9F7A68", 
               "#4E2901", "#A9A7AE", "#4C4A4B", "#000000")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})
kp.stro <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.stro[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#8B0000", "#CA4B2C", "#F8766D", "#DE8C00", 
               "#D8C920", "#A8C709", "#094500", "#00BA38", 
               "#1DDED0", "#224AD9", "#000168", "#6E15D4", 
               "#C77CFF", "#FB61DE", "#C0BC87", "#9F7A68", 
               "#4E2901", "#A9A7AE", "#4C4A4B", "#000000")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})
kp.neur <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.neur[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#8B0000", "#CA4B2C", "#F8766D", "#DE8C00", 
               "#D8C920", "#A8C709", "#094500", "#00BA38", 
               "#1DDED0", "#224AD9", "#000168", "#6E15D4", 
               "#C77CFF", "#FB61DE", "#C0BC87", "#9F7A68", 
               "#4E2901", "#A9A7AE", "#4C4A4B", "#000000")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})
kp.trop <- lapply(sen_genes, function(g) {
  p <- plot_genes_in_pseudotime(
    cds_subset = cds.subset.trop[g,], 
    color_cells_by = "Timepoint"
  ) + theme(
    text = element_text(size = 30), 
    legend.key.size = unit(0.75, "cm")
  ) + scale_color_manual(
    values = c("#8B0000", "#CA4B2C", "#F8766D", "#DE8C00", 
               "#D8C920", "#A8C709", "#094500", "#00BA38", 
               "#1DDED0", "#224AD9", "#000168", "#6E15D4", 
               "#C77CFF", "#FB61DE", "#C0BC87", "#9F7A68", 
               "#4E2901", "#A9A7AE", "#4C4A4B", "#000000")
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
    values = c("#8B0000", "#CA4B2C", "#F8766D", "#DE8C00", 
               "#D8C920", "#A8C709", "#094500", "#00BA38", 
               "#1DDED0", "#224AD9", "#000168", "#6E15D4", 
               "#C77CFF", "#FB61DE", "#C0BC87", "#9F7A68", 
               "#4E2901", "#A9A7AE", "#4C4A4B", "#000000")
  )
  p + coord_cartesian(ylim = ylims[[g]])
})

png(
  filename = "../../Images/MouseReprogramming/scbg_kineticplot_somatic.png", 
  width = 948, height = 1000, res = 120
)
wrap_plots(
  kp.som, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_kineticplot_stromal.png", 
  width = 948, height = 1000, res = 120
)
wrap_plots(
  kp.stro, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_kineticplot_neural.png", 
  width = 940, height = 1000, res = 120
)
wrap_plots(
  kp.neur, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_kineticplot_trophoblast.png", 
  width = 930, height = 1000, res = 120
)
wrap_plots(
  kp.trop, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_kineticplot_pluripotent.png", 
  width = 1100, height = 1000, res = 120
)
wrap_plots(
  kp.plu, ncol = 1, guides = "collect", 
  axes = "collect", axis_titles = "collect"
)
dev.off()


