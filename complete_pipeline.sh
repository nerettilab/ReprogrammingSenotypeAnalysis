#!/bin/bash

# Initialize conda (must add location of conda directory)
source  /<your_directory_here>/miniforge3/etc/profile.d/conda.sh

### Conda environment must be created before rest of the code will work.
### If you have already created this environment, you can skip this part.
conda env create -f CondaEnvironment/init_analysis_env.yml
conda deactivate
conda activate analysis_env
Rscript CondaEnvironment/init_analysis_env.R
###
###

# Activate conda environment
conda deactivate
conda activate analysis_env

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

# Download data from publically available data sets of senescence induction
humandir="Data/SenotypeComparisons"

wget -O $humandir/BJ_SIS_alspach2014_young_9_11.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358449/suppl/GSM1358449_CCGATTA_gene_counts.txt.gz
wget -O $humandir/BJ_SIS_alspach2014_young_9_24.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358452/suppl/GSM1358452_TACTCTA_gene_counts.txt.gz
wget -O $humandir/BJ_SIS_alspach2014_sen_9_11.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358451/suppl/GSM1358451_GGCAGCG_gene_counts.txt.gz
wget -O $humandir/BJ_SIS_alspach2014_sen_9_24.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358447/suppl/GSM1358447_AGACTGA_gene_counts.txt.gz

wget -O $humandir/BJ_OIS_loayza-puch2013_FPKM.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42509/suppl/GSE42509_PQST_rna_rp_fpkms.txt.gz

wget -O $humandir/BJ_RS_marthandam2016_RPKM.xls.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63577/suppl/GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz

gunzip $humandir/*.gz

# Collapse single-cell reprogramming data sets into pseudo-bulk for comparisons with senescence data
Rscript Scripts/SenotypeComparisons/collapse_into_pseudobulk.R

# Compare senotypes of reprogramming data with those of senescence induction data
Rscript Scripts/SenotypeComparisons/senescence_gene_comparisons.R
Rscript Scripts/SenotypeComparisons/senotype_rank-correlation.R
