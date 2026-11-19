library(remotes)

# Update to correct version of glmGamPoi. Needed for normalization via SCTransform
install_version("glmGamPoi", version = "1.20.0", upgrade = "never", repos = "https://bioconductor.org/packages/3.21/bioc")

# Needed for identifying conserved cell type markers
install_version("metap", version = "1.12", upgrade = "never", repos = "https://cloud.r-project.org")

# Needed for finding markers of RNA clusters
remotes::install_github("immunogenomics/presto", ref = "7636b3d", upgrade = "never")

# Needed for analyzing ECM enrichment in data
remotes::install_github("Matrisome/MatrisomeAnalyzeR", ref = "d8aec6d", upgrade = "never")

# Install DoubletFinder
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", ref = "1b244d8", upgrade = "never")

# For retrieving datasets from Seurat vignettes
remotes::install_github("satijalab/seurat-data", ref = "3e51f44", upgrade = "never")

# Install SeuratWrappers for additional functionality between Seurat objects and Monocle
remotes::install_github("satijalab/seurat-wrappers", ref = "050f4f0", upgrade = "never")
