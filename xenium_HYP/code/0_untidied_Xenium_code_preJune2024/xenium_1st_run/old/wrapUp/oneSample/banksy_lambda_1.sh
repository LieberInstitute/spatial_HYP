#!/bin/sh

#SBATCH --partition=shared 

#SBATCH --cpus-per-task=1

ml conda_R/4.3

R CMD BATCH banksy_lambda_1.R



