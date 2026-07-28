# Reprogramming Senotype Analysis

Nunez and Neretti (in preparation)

A distinct senotype emerges during OSKM-mediated cellular reprogramming

This repository contains the raw data and code required to set up Conda Environment and create the figures used in this study. Destination folders will be created for the data that will be generated.

To install and run:
1. Clone the GitHub repository
2. Open CondaEnvironment/init_analysis_env.yml and edit the last line to input location of your Conda directory
3. Open complete_pipeline.sh and edit line 4 to input location of your Conda directory
4. Run complete_pipeline.sh

Running complete_pipeline.sh should automatically create/activate the Conda Environment and run all of the necessary scripts. If complete_pipeline.sh does not function properly, then the following steps may be followed instead to achieve the intended result:

1. Initialize Conda
2. Create Conda environment with: conda env create -f CondaEnvironment/init_analysis_env.yml
3. Activate Conda envrionment, then download remaining dependencies with: Rscript CondaEnvironment/init_analysis_env.R
4. Run the following scripts in the following order:
   - ./Scripts/download_data.sh
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
   - Rscript Scripts/SenotypeComparisons/collapse_into_pseudobulk.R
   - Rscript Scripts/SenotypeComparisons/senescence_gene_comparisons.R