#!/bin/bash
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH -J msHYP_to_SCE
#SBATCH -o 01-jhpce_logs/mshyp_to_SCE.out
#SBATCH -e 01-jhpce_logs/mshyp_to_SCE.err

module load conda_R/4.3.x

# dispatch from directory /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/analysis/06-BayesSpace

R CMD BATCH --no-save --no-restore 01a-convert_Yao23HYP_to_SCE.R