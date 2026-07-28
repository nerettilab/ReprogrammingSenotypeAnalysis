#!/bin/bash

# Initialize conda (must add location of conda directory)
source  /<your_directory_here>/miniforge3/etc/profile.d/conda.sh


## Conda environment must be created before rest of the code will work.
## If you have already created this environment, you can skip this part.
./CondaEnvironment/init_analysis_env.sh


# Activate conda environment
conda deactivate
conda activate analysis_env

# Download necessary data and create destination directories
./Scripts/download_data.sh

# Create Seurat objects for human reprogramming data sets
Rscript Scripts/Xing2020/create_seurat_objects.R
Rscript Scripts/Panariello2023/create_seurat_objects.R

# Create Seurat objects for mouse reprogramming data sets
Rscript Scripts/Schiebinger2019/create_seurat_objects.R

# Merge human reprogramming Seurat objects
Rscript Scripts/MergedHuman/merge_seurat_objects.R

# Monocle3 trajectory analysis of human and mouse reprogramming data
Rscript Scripts/MergedHuman/monocle3_analysis.R
Rscript Scripts/Schiebinger2019/monocle3_analysis.R

# Compare expression of cell fate markers between trajectory markers
Rscript Scripts/MergedHuman/trajectory_markers.R
Rscript Scripts/Schiebinger2019/trajectory_markers.R

# Analyze expression of senescence genes in reprogramming data
Rscript Scripts/MergedHuman/senescence_gene_analysis.R
Rscript Scripts/Schiebinger2019/senescence_gene_analysis.R

# Identify senescence markers shared between human and mouse reprogramming data
Rscript Scripts/HumanMouseComparisons/find_common_senescence_markers.R

# Collapse single-cell reprogramming data sets into pseudo-bulk for comparisons with senescence data
Rscript Scripts/SenotypeComparisons/collapse_into_pseudobulk.R

# Compare senotypes of reprogramming data with those of senescence induction data
Rscript Scripts/SenotypeComparisons/senescence_gene_comparisons.R
