#!/bin/bash
#SBATCH --mem=10G
#SBATCH --job-name=04samui
#SBATCH --output=logs/%x%a.txt
#SBATCH --array=1-13

# -4%4

# #SBATCH --error=logs/%x%a.txt

#echo -e "\n"
echo "date: $(date)"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load samui/1.0.0-next.49
python 03_loopy-xenium.py

echo "**** Job ends ****"
date

#echo -e "\n"
