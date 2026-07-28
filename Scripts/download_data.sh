#!/bin/bash

## Download and counts data and metadata from Xing 2020
# Create destination directory
mkdir -p Data/Xing2020/GSE118258
xingdir="Data/Xing2020/GSE118258"

# Download and decompress data
wget -O $xingdir/GSE118258_Annotation.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118258/suppl/GSE118258_Annotation.txt.gz
wget -O $xingdir/GSE118258_UMI.csv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118258/suppl/GSE118258_UMI.csv.gz
gunzip $xingdir/*.gz


## Download and counts data and metadata from Panariello 2023
# Create destination directories
mkdir -p Data/Panariello2023/GSE221739 Data/Panariello2023/seurat_objects_raw Data/Panariello2023/seurat_objects_filtered
panadir="Data/Panariello2023/GSE221739/"
mkdir $panadir/D00 $panadir/D03 $panadir/D05 $panadir/D07 $panadir/D09 $panadir/D11 $panadir/D13 $panadir/D15

# Download and decompress data
wget -O $panadir/GSE221739_RAW.tar https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221739/suppl/GSE221739_RAW.tar
tar -C $panadir/D00 -xf $panadir/GSE221739_RAW.tar GSM6894025_* GSM6894026_*
tar -C $panadir/D03 -xf $panadir/GSE221739_RAW.tar GSM6894027_* GSM6894028_*
tar -C $panadir/D05 -xf $panadir/GSE221739_RAW.tar GSM6894029_* GSM6894032_*
tar -C $panadir/D07 -xf $panadir/GSE221739_RAW.tar GSM6894033_* GSM6894035_*
tar -C $panadir/D09 -xf $panadir/GSE221739_RAW.tar GSM6894036_* GSM6894039_*
tar -C $panadir/D11 -xf $panadir/GSE221739_RAW.tar GSM6894040_* GSM6894042_*
tar -C $panadir/D13 -xf $panadir/GSE221739_RAW.tar GSM6894043_* GSM6894044_*
tar -C $panadir/D15 -xf $panadir/GSE221739_RAW.tar GSM6894045_* GSM6894046_*
rm $panadir/GSE221739_RAW.tar
gunzip $panadir/*/*


## Create destination directory to be used when human data sets are merged
mkdir Data/MergedHuman


## Download and counts data and metadata from Schiebinger 2019
# Create destination directory
mkdir -p Data/Schiebinger2019/GSE115943
scbgdir="Data/Schiebinger2019/GSE115943"

# Download and decompress data
wget -O $scbgdir/GSE115943_RAW.tar https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115943/suppl/GSE115943_RAW.tar
tar -C $scbgdir -xzf $scbgdir/GSE115943_RAW.tar
rm GSE115943_RAW.tar


## Download data from publically available data sets of senescence induction
# Create destination directory if it does not exist
sencomdir="Data/SenotypeComparisons"

# Download data
# Stress-induced senescence
wget -O $sencomdir/BJ_SIS_alspach2014_young_9_11.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358449/suppl/GSM1358449_CCGATTA_gene_counts.txt.gz
wget -O $sencomdir/BJ_SIS_alspach2014_young_9_24.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358452/suppl/GSM1358452_TACTCTA_gene_counts.txt.gz
wget -O $sencomdir/BJ_SIS_alspach2014_sen_9_11.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358451/suppl/GSM1358451_GGCAGCG_gene_counts.txt.gz
wget -O $sencomdir/BJ_SIS_alspach2014_sen_9_24.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358447/suppl/GSM1358447_AGACTGA_gene_counts.txt.gz
# Oncogene-induced senescence
wget -O $sencomdir/BJ_OIS_loayza-puch2013_FPKM.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42509/suppl/GSE42509_PQST_rna_rp_fpkms.txt.gz
# Replicative senescence
wget -O $sencomdir/BJ_RS_marthandam2016_RPKM.xls.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63577/suppl/GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz

gunzip $sencomdir/*.gz # Decompress data

## SIDE NOTE about GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz:
## This Excel spreadsheet has two pages - one for raw counts, and 
## one for RPKM; only the "counts" page will be needed. This page 
## has been isolated an converted into a CSV file for ease of use 
## in R. The CSV is included in the repository.


# Make directories for images that will be generated
mkdir -p Images/HumanReprogramming Images/MouseReprogramming Images/HumanMouseComparisons Images/SenotypeComparisons