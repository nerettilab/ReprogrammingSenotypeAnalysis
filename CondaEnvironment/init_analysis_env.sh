#!/bin/bash
#SBATCH -n 16
#SBATCH --mem=200G
#SBATCH -t 48:00:00
#SBATCH -o init_env-output-%j.out
#SBATCH -e init_env-output-%j.err

# Initialize conda
source /<your_directory_here>/miniforge3/etc/profile.d/conda.sh

# Create conda environment from yml file
conda env create -f init_analysis_env.yml
conda deactivate
conda activate analysis_env

# Run RScript to install other dependencies for environment in R
Rscript init_analysis_env.R
