#!/bin/bash
#SBATCH --mem=50G
#SBATCH --job-name=r2py-xen
#SBATCH --array=1

echo -en '\n'
echo "**** Job starts ****"
date
echo -en '\n'

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.4
Rscript 01_r2py-xen.R

mkdir -p ./logs/

cat ./slurm-${SLURM_JOBID}*.out >> ./logs/01_r2py-xe.txt
echo -en '\n'
rm ./slurm-${SLURM_JOBID}*.out

echo "**** Job ends ****"
date