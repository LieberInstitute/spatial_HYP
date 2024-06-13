#!/bin/sh

#SBATCH --partition=shared 

#SBATCH --cpus-per-task=1

ml conda_R/4.3

R CMD BATCH step4_findMarkerGenes_BanksyClusters_lambda0.R





