#!/bin/bash
#SBATCH --mem=30G
#SBATCH -n 1
#SBATCH --job-name=hyp_bsk2
#SBATCH -o jhpce_logs/hyp_bsk2.txt
#SBATCH -e jhpce_logs/hyp_bsk2.txt
#SBATCH --array=1-12

module load conda_R/4.3.x

R CMD BATCH 01-BayesSpace_hvg1020_nnsvg1020_harmony_k2.R
