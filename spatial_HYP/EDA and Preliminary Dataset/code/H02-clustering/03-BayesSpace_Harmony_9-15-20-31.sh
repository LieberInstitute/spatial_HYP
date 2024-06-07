#
#$ -l h_fsize=20G,mem_free=25G,h_vmem=38G 
#$ -cwd
#$ -m e
#$ -t 1-16

module load conda_R/4.2

R CMD BATCH 03-BayesSpace_hvg1020_nnsvg1020_harmony_k_9-15-20-31.R
