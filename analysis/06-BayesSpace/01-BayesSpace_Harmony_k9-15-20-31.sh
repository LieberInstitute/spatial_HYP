#!/bin/bash
#SBATCH --mem=45G
#SBATCH -n 1
#SBATCH --job-name=hyp_bspaces
#SBATCH -o 01-jhpce_logs/hyp_bsk9-15-20-31.out
#SBATCH -e 01-jhpce_logs/hyp_bsk9-15-20-31.err
#SBATCH --array=1-48

# run from the directory /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/analysis/06-BayesSpace

module load conda_R/4.3.x

R CMD BATCH 01-BayesSpace_hvg1020_nnsvg1020_harmony_k9-15-20-31.R
