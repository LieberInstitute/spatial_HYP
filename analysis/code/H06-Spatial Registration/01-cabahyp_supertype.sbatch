#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH -J abahyp_subcl
#SBATCH -e logs/err_supertype.out
#SBATCH -o logs/log_supertype.out

cd /users/bmulvey/spatreg_hyp_ABAhyp

module load R/4.3

R CMD BATCH --no-save --no-restore reg_wrapper_ABAhyp_supertype.R
