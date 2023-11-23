#!/bin/bash
#SBATCH --mem=45G
#SBATCH -n 1
#SBATCH --job-name=hyp_bsk2
#SBATCH -o jhpce_logs/hyp_bsk9-15-20-31.txt
#SBATCH -e jhpce_logs/hyp_bsk9-15-20-31.txt
#SBATCH --array=1-48

module load conda_R/4.3.x

R CMD BATCH 01-BayesSpace_hvg1020_nnsvg1020_harmony_k9-15-20-31.R
