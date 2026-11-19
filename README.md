# Reprogramming Senotype Analysis

Nunez and Neretti (Accepted - NPJ Aging)

A distinct senotype emerges during OSKM-mediated cellular reprogramming

This repository contains the raw data and code required to set up Conda Environment and create the figures used in this study. Destination folders will be created for the data that will be generated.

To install and run:
1. Clone the GitHub repository
2. Open CondaEnvironment/init_analysis_env.yml and edit the last line to input location of your Conda directory
3. Open complete_pipeline.sh and edit line 4 to input location of your Conda directory
4. Run complete_pipeline.sh

Running complete_pipeline.sh should automatically create/activate the Conda Environment and run all of the necessary scripts. If complete_pipeline.sh does not function properly, then the following steps may be used instead to achieve the intended result:

1. Initialize Conda
2. Create Conda environment with: conda env create -f CondaEnvironment/init_analysis_env.yml
3. Activate Conda envrionment, then download remaining dependencies with: Rscript CondaEnvironment/init_analysis_env.R
4. Run the following scripts in the following order:
   - Rscript Scripts/Xing2020/create_seurat_objects.R
   - Rscript Scripts/Panariello2023/create_seurat_objects.R
   - Rscript Scripts/Schiebinger2019/create_seurat_objects.R
   - Rscript Scripts/MergedHuman/merge_seurat_objects.R
   - Rscript Scripts/MergedHuman/monocle3_analysis.R
   - Rscript Scripts/Schiebinger2019/monocle3_analysis.R
   - Rscript Scripts/MergedHuman/trajectory_markers.R
   - Rscript Scripts/Schiebinger2019/trajectory_markers.R
   - Rscript Scripts/MergedHuman/senescence_gene_analysis.R
   - Rscript Scripts/Schiebinger2019/senescence_gene_analysis.R
   - Rscript Scripts/HumanMouseComparisons/find_common_senescence_markers.R
5. Run the following lines of code:
   wget -O Data/SenotypeComparisons/BJ_SIS_alspach2014_young_9_11.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358449/suppl/GSM1358449_CCGATTA_gene_counts.txt.gz
   wget -O Data/SenotypeComparisons/BJ_SIS_alspach2014_young_9_24.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358452/suppl/GSM1358452_TACTCTA_gene_counts.txt.gz
   wget -O Data/SenotypeComparisons/BJ_SIS_alspach2014_sen_9_11.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358451/suppl/GSM1358451_GGCAGCG_gene_counts.txt.gz
   wget -O Data/SenotypeComparisons/BJ_SIS_alspach2014_sen_9_24.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1358nnn/GSM1358447/suppl/GSM1358447_AGACTGA_gene_counts.txt.gz
   wget -O Data/SenotypeComparisons/BJ_OIS_loayza-puch2013_FPKM.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42509/suppl/GSE42509_PQST_rna_rp_fpkms.txt.gz
   wget -O Data/SenotypeComparisons/BJ_RS_marthandam2016_RPKM.xls.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63577/suppl/GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz
   gunzip Data/SenotypeComparisons/*.gz
6. Run the following scripts in the following order:
   - Rscript Scripts/SenotypeComparisons/collapse_into_pseudobulk.R
   - Rscript Scripts/SenotypeComparisons/senescence_gene_comparisons.R
   - Rscript Scripts/SenotypeComparisons/senotype_rank-correlation.R
