#!/bin/bash
#SBATCH --mem=285G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -o ./logs/xeniumresegment_2023slides_%a.out
#SBATCH -e ./logs/xeniumresegment_2023slides_%a.err
#SBATCH --array=1-6
#SBATCH --constraint="intel"

echo "**** Job starts ****"
date

# must be run from code/01_
# and files must be moved from code/
# 01_ after xeniumranger finishes
# because XR won't allow slashes in the output directory name (--id) argument

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module --ignore_cache load xeniumranger/1.7.1

# get current sample to resegment
FILE=$( sed -n ${SLURM_ARRAY_TASK_ID}p xeniumranger_sampledirs.txt)

# predefined output dirs
#  which will be in the code directory temporarily
OUTPATH=$( sed -n ${SLURM_ARRAY_TASK_ID}p reseg_outdirs.txt)

# get full path so we can relocate output without problems if
# movedd to different dcs0N later
FULLPATH=$(realpath ../..)

# run xeniumranger, with localcores set to 
# the number of cores in the sbatch job
# and, since using one node, localmem set to 90% 
# of the SBATCH requested memory
# both per xeniumranger documentation. (otherwise,
#  xenranger will expect it has 1TB of
#  RAM at its disposal and crash when it ineveitably doesnt)
xeniumranger resegment --xenium-bundle=${FILE} --resegment-nuclei=false --id=${OUTPATH} --localcores=48 --localmem=250 --disable-ui=true --jobmode=local

# move output where it should've been able to go in the first place
mv ${OUTPATH} ${FULLPATH}/processed-data/01_xeniumranger1.7-resegment/
echo "**** Job ends ****"
date
