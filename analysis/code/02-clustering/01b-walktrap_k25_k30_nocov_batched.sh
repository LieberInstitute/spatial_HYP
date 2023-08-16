#
#$ -l h_fsize=20G,mem_free=25G,h_vmem=38G 
#$ -cwd
#$ -m e
#$ -t 1-14

module load conda_R/4.2

R CMD BATCH walktrap_k2530_052623_batched.R
