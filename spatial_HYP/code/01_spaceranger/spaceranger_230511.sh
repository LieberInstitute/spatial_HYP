#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 8
#$ -N spatialHYP_spaceranger_230511
#$ -o logs/spaceranger_230511.$TASK_ID.txt
#$ -e logs/spaceranger_230511.$TASK_ID.txt
#$ -m e
#$ -t 1-8
#$ -tc 8

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load SpaceRanger
module load spaceranger/2.0.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" samples_230511.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
SAM=$(paste <(echo ${SLIDE}) <(echo "-") <(echo ${CAPTUREAREA}) -d '')
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH=$(ls -d /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/raw-data/FASTQ/${SAMPLE}/)

## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH}\
    --image=/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/processed-data/Images/${SAMPLE}.tif \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --loupe-alignment=/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/processed-data/Images/loupe/${SAM}.json \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving results to new location"
date
mkdir -p /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/processed-data/01_spaceranger/spaceranger_230511/
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/processed-data/01_spaceranger/spaceranger_230511/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
