#!/bin/bash
#SBATCH --job-name=BNKS_HYP_mult
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=250G
#SBATCH --constraint="znver2|znver3|intel"
#SBATCH -o ./logs/01-BanksyM0lam0_leidenRes2_multi.out
#SBATCH -e ./logs/01-BanksyM0lam0_leidenRes2_multi.err


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3.x

## run this script from xenium_HYP/code/05blahblah
Rscript 01-Banksy_M0lam0_leidenres2_multisamp.R

echo "**** Job ends ****"
date

