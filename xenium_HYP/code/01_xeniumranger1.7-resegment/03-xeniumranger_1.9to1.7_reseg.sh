#!/bin/bash
#SBATCH --mem=285G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH -o ./logs/xeniumresegment_2024slides_%a.out
#SBATCH -e ./logs/xeniumresegment_2024slides_%a.err
#SBATCH --array=1-7
#SBATCH --constraint="(amd&znver2)|(amd&znver3)|intel"

echo "**** Job starts ****"
date

# required argument for AMD nodes to proceed w xenranger;
# default (soft) limit is 1024 on such nodes.
# see https://www.10xgenomics.com/support/software/xenium-ranger/latest/release-notes/XR-system-requirements#resource-limits
prlimit --nofile=16000

# must be run from code/01_
# and files must be moved from code/01_ after xeniumranger finishes
# because XR won't allow slashes in the output directory name (--id) argument

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module --ignore_cache load xeniumranger/1.7.1

# get current sample to resegment
FILE=$( sed -n ${SLURM_ARRAY_TASK_ID}p 03-xenranger1.9_sampledirs.txt)

# predefined output dirs
#  which will be in the code directory temporarily
OUTPATH=$( sed -n ${SLURM_ARRAY_TASK_ID}p 03-reseg_outdirs.txt)

# get full path so we can relocate output without problems if
# movedd to different dcs0N later
FULLPATH=$(realpath ../..)

# run xeniumranger, with localcores set to 
# the number of cores in the sbatch job *2 (intel) or *1 (amd)
# and, since using one node, localmem set to 90% 
# of the SBATCH requested memory
# both per xeniumranger documentation. (otherwise,
#  xenranger will expect it has 1TB of
#  RAM at its disposal and crash when it ineveitably doesnt)
xeniumranger resegment --xenium-bundle=${FILE} --resegment-nuclei=true --id=${OUTPATH} --localcores=36 --localmem=255 --disable-ui=true --jobmode=local

# move output where it should've been able to go in the first place
mv ${OUTPATH} ${FULLPATH}/processed-data/01_xeniumranger1.7-resegment/
echo "**** Job ends ****"
date
