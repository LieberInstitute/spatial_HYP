#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task=6
#SBATCH -J abahyp_subcl
#SBATCH -e logs/err_subclass.out
#SBATCH -o logs/log_subclass.out

cd /users/bmulvey/spatreg_hyp_ABAhyp

module load R/4.3

R CMD BATCH --no-save --no-restore reg_wrapper_ABAhyp_subclass.R
