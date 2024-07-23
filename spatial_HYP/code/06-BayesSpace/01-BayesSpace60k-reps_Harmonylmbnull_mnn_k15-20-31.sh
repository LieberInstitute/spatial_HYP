#!/bin/bash
#SBATCH --mem=136G
#SBATCH --cpus-per-task=2
#SBATCH -J hyp_bspaces
#SBATCH -o 01-jhpce_logs/hyp_bs50kiter_k-15-20-31.out
#SBATCH -e 01-jhpce_logs/hyp_bs50kiter_k-15-20-31.err
#SBATCH --array=1-18

module load conda_R/4.3.x

# dispatch from directory /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/spatial_HYP/code/06-BayesSpace

R CMD BATCH --no-save --no-restore 01-BayesSpace-10k-burn-50k-addnl-reps_hvg20_nnsvg10nnsvg20_harmony_k15-20-31.R

