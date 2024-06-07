#!/bin/bash
#SBATCH --cpus-per-task 4
#SBATCH --mem=40G
#SBATCH -J make_HYPspe
#SBATCH -o logs/makespe_wV13M13-362.out #standard out goes here
#SBATCH -e logs/makespe_wV13M13-362.err #standard error goes here

echo "**** Job starts ****"
date

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02b_raw_spe_w-V13M13-362.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
