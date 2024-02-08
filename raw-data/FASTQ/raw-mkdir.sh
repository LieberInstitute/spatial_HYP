#!/bin/bash
#SBATCH --array=1

SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" raw-mkdir.txt)

mkdir -p ./logs/
mkdir -p /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/raw-data/FASTQ/${SAMPLE}/

mv slurm-*.out ./logs