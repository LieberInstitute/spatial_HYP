#!/bin/sh

#SBATCH --partition=shared 

#SBATCH --cpus-per-task=1

ml conda_R/4.3

R CMD BATCH joint_runs_multisample_Banksy_lambda0_res2.R





