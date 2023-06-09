#
#$ -l h_fsize=10G,mem_free=10G,h_vmem=25G
#$ -pe local 8
#$ -cwd
#$ -m e
#$ -t 1-16

cd /users/bmulvey/precast-hvgsvg

module load conda_R/4.2

R CMD BATCH 04b-precast_hvg-svg_9-15-20-31.R
