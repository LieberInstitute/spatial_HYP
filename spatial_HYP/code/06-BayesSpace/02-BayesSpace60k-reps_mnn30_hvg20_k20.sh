#!/bin/bash
#SBATCH --mem=125G
#SBATCH --cpus-per-task=4
#SBATCH --constraint="intel"
#SBATCH -t 1-12:00
#SBATCH -J hypBS20_mnnhvg
#SBATCH -o 02-jhpce_logs/hyp_bs50kiter_mnn30_hvg20_k20.out
#SBATCH -e 02-jhpce_logs/hyp_bs50kiter_mnn30_hvg20_k20.err

module load conda_R/4.3.x

# dispatch from directory /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/spatial_HYP/code/06-BayesSpace

R CMD BATCH --no-save --no-restore 02-BayesSpace-10k-burn-50k-addnl-reps_hvg20_mnn30_k20.R

