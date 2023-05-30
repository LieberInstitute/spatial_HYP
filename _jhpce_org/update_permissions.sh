#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=400G
#$ -N update_permissions
#$ -o logs/update_permissions.txt
#$ -e logs/update_permissions.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


sh /dcs04/lieber/lcolladotor/_jhpce_org_LIBD001/update_permissions_spatialteam.sh /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/
