######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

lambda=0
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed_banksy =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"

###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out_compare_cell_nucleus/"
dir.create(dir_out)
setwd(dir_out)

###########
dir_out_nu = "/dcs04/hansen/data/ywang/xenium_joint_runs/out_nucleus/"
dir_out_cell = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"


###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
aname <- "normcounts"

###########
sample_names_old =  c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")
sample_names_new =   c("br1225a","br1225b","br8741c",
                       "br8741d","br5459a","br5459b",
                       "br8667c")

###########
res=2
se_cell = readRDS(paste0(dir_out_cell, "/spe_joint_2runs_banksy_lambda0_res", res, ".rds"))
se_nucleus = readRDS(paste0(dir_out_nu, "/spe_joint_2runs_banksy_lambda0_res", res, ".rds"))
dim(se_nucleus)
dim(se_cell)
# [1]    366 879396
# [1]    366 879396

cellClu = colData(se_cell)[,"clust_M0_lam0_k50_res2"]

############
# counts_nu = assays(se_nucleus)$counts
sampleID_nu = colData(se_nucleus)[,"sample_id"]
cellClu = colData(se_cell)[,"clust_M0_lam0_k50_res2"]

############
total = c()
total_ESR1expressingOly = c()
total_ESR1TAC3expressingOly = c()

both_GAD1_SLC17A6 = c()
both_GAD1_SLC17A6_ESR1expressingOly = c()
both_GAD1_SLC17A6_ESR1TAC3expressingOly = c()



gad1Only = c()
gad1Only_ESR1expressingOly = c()
gad1Only_ESR1TAC3expressingOly = c()

SLC17A6only = c()
SLC17A6only_ESR1expressingOly = c()
SLC17A6only_ESR1TAC3expressingOly = c()

sample_unique = unique(sampleID_nu)

for(sample_tmp in sample_unique){
  print(sample_tmp)
  which_cells_tmp = which(sampleID_nu == sample_tmp & cellClu=="33")
  counts_nu_tmp = assays(se_nucleus)$counts[,which_cells_tmp]
  
  ESR1 = counts_nu_tmp["ESR1",]
  GAD1 = counts_nu_tmp["GAD1",]
  SLC17A6 = counts_nu_tmp["SLC17A6",]
  TAC3 = counts_nu_tmp["TAC3",]
  
  both_GAD1_SLC17A6 = c(both_GAD1_SLC17A6,sum(GAD1>0 & SLC17A6>0))  
  gad1Only = c(gad1Only,sum(GAD1>0 & SLC17A6==0))  
  SLC17A6only = c(SLC17A6only,sum(GAD1==0 & SLC17A6>0))  
  total = c(total, ncol(counts_nu_tmp))
  
  both_GAD1_SLC17A6_ESR1expressingOly = c(both_GAD1_SLC17A6_ESR1expressingOly,sum(ESR1>0 & (GAD1>0 & SLC17A6>0))  )
  gad1Only_ESR1expressingOly = c(gad1Only_ESR1expressingOly,sum(ESR1>0 & (GAD1>0 & SLC17A6==0))  )
  SLC17A6only_ESR1expressingOly = c(SLC17A6only_ESR1expressingOly,sum( ESR1>0 & (GAD1==0 & SLC17A6>0)) ) 
  total_ESR1expressingOly = c(total_ESR1expressingOly, sum(ESR1>0))
  
  both_GAD1_SLC17A6_ESR1TAC3expressingOly = c(both_GAD1_SLC17A6_ESR1TAC3expressingOly,sum((ESR1>0&TAC3>0) & (GAD1>0 & SLC17A6>0))  )
  gad1Only_ESR1TAC3expressingOly = c(gad1Only_ESR1TAC3expressingOly,sum((ESR1>0&TAC3>0) & (GAD1>0 & SLC17A6==0))  )
  SLC17A6only_ESR1TAC3expressingOly = c(SLC17A6only_ESR1TAC3expressingOly,sum( (ESR1>0&TAC3>0) & (GAD1==0 & SLC17A6>0)) ) 
  total_ESR1TAC3expressingOly = c(total_ESR1TAC3expressingOly, sum(ESR1>0&TAC3>0))
}
# sample_unique
# [1] "br6197"  "br5993a" "br5993b" "br1735a" "br1735b" "br6588"  "br1225a"
# [8] "br1225b" "br8741c" "br8741d" "br5459a" "br5459b" "br8667c"

tb_ = rbind(both_GAD1_SLC17A6, gad1Only, SLC17A6only, total)
tb_esr1 = rbind(both_GAD1_SLC17A6_ESR1expressingOly, gad1Only_ESR1expressingOly, SLC17A6only_ESR1expressingOly, total_ESR1expressingOly)
tb_esr1_tac3 = rbind(both_GAD1_SLC17A6_ESR1TAC3expressingOly, gad1Only_ESR1TAC3expressingOly, SLC17A6only_ESR1TAC3expressingOly, total_ESR1TAC3expressingOly)
rownames(tb_) =rownames(tb_esr1) =rownames(tb_esr1_tac3) = c("both_GAD1_SLC17A6", "gad1Only", "SLC17A6only"," total")
colnames(tb_) =colnames(tb_esr1) =colnames(tb_esr1_tac3) = sample_unique
write.csv(tb_, file = "tb_GAD1_SLC17A6_clu33.csv")
write.csv(tb_esr1, file = "tb_GAD1_SLC17A6_clu33_ESR1expressingOnly.csv")
write.csv(tb_esr1_tac3, file = "tb_GAD1_SLC17A6_clu33_ESR1TAC3expressingOnly.csv")

# ###########
# spe_joint_nu_old = readRDS(paste0(dir_data_processed_downstream_old_nucleud, "xens.rds"))
# spe_joint_cell_old = readRDS(paste0(dir_data_processed_downstream_old, "xens.rds"))
# 
# spe_joint_nu_new = readRDS(paste0(dir_data_processed_downstream_new_nucleud, "xens.rds"))
# spe_joint_cell_new = readRDS(paste0(dir_data_processed_downstream_new, "xens.rds"))
# 

head(colData(spe_joint_nu_old))
head(colData(spe_joint_nu_new))
head(colData(spe_joint_cell_old))
head(colData(spe_joint_cell_new))

###########
colData(spe_joint_nu_old)
# colData(spe_joint_cell)

counts_cell_old =assays(spe_joint_cell_old)$counts
counts_cell_new =assays(spe_joint_cell_new)$counts

counts_nu_old =assays(spe_joint_nu_old)$counts
counts_nu_new =assays(spe_joint_nu_new)$counts

colnames(counts_cell_old) = sub(".*[.]","", colnames(counts_cell_old))
colnames(counts_cell_new) = sub(".*[.]","", colnames(counts_cell_new))
colnames(counts_nu_old) = sub(".*[.]","", colnames(counts_nu_old))
colnames(counts_nu_new) = sub(".*[.]","", colnames(counts_nu_new))


length(intersect(colnames(counts_cell_old), colnames(counts_nu_old)))
length(intersect(colnames(counts_cell_new), colnames(counts_nu_new)))
# length(intersect(colnames(counts_nu_old), colnames(counts_cell_new)))
colnames(counts_nu_new)[1]
colnames(counts_cell_new)[1]

#########

shared_genes_cell = intersect(rownames(counts_cell_old),
                         rownames(counts_cell_new))

shared_genes_nu = intersect(rownames(counts_nu_old),
                         rownames(counts_nu_new))

shared_genes = intersect(shared_genes_cell, shared_genes_nu)

counts_cell = cbind(counts_cell_old[shared_genes, ],
                    counts_cell_new[shared_genes,])
counts_nu  = cbind(counts_nu_old[shared_genes, ],
                   counts_nu_new[shared_genes,])
dim(counts_nu)
# > length(cells_shared)
# [1] 408338
# > dim(counts_nu)
# [1]    434 838037

###########
# 

counts_ESR1_cell = counts_cell["ESR1",]
counts_ESR1_nu = counts_nu["ESR1",]

genes_shared  =  intersect(rownames(counts_cell),
                            rownames(counts_nu))

cells_shared = intersect(colnames(counts_cell),
                         colnames(counts_nu))

length(cells_shared)

set.seed(111)
cells_shared_subsample = sample(cells_shared, 1000)
# pdf("")

####### make comparison plots
library(scales)

# plot(counts_ESR1_cell[cells_shared_subsample],col = alpha("black",.01), counts_ESR1_nu[cells_shared_subsample], xlab = "cell", ylab="nuclear", main = "ESR1 gene counts per cell, color transparency = 0.01")
# abline(0,1,col = alpha("blue",.5))

plot(counts_ESR1_cell[cells_shared_subsample],col = alpha("black",.1), counts_ESR1_nu[cells_shared_subsample], xlab = "cell", ylab="nuclear", main = "ESR1 gene counts per cell, color transparency = 0.1")
abline(0,1,col = alpha("blue",.5))

# plot(counts_ESR1_cell[cells_shared_subsample],col = alpha("black",1), counts_ESR1_nu[cells_shared_subsample], xlab = "cell", ylab="nuclear", main = "ESR1 gene counts per cell, color transparency = 1")
# abline(0,1,col = alpha("blue",.5))

####### total counts
csum_cell = colSums(counts_cell[genes_shared,])
csum_nu = colSums(counts_nu[genes_shared,])

# plot(csum_cell[cells_shared_subsample],col = alpha("black",.01), csum_nu[cells_shared_subsample], xlab = "cell", ylab="nuclear", main = "total gene counts per cell, color transparency = 0.01")
# abline(0,1,col = alpha("blue",.5))
# 
plot(csum_cell[cells_shared_subsample],col = alpha("black",.1), csum_nu[cells_shared_subsample], xlab = "cell", ylab="nuclear", main = "total gene counts per cell, color transparency = 0.1")
abline(0,1,col = alpha("blue",.5))

sum(csum_cell[cells_shared])
sum(csum_nu[cells_shared])

# > sum(csum_cell[cells_shared])
# [1] 308845512
# > sum(csum_nu[cells_shared])
# [1] 80851215

sum(csum_cell[cells_shared_subsample])
sum(csum_nu[cells_shared_subsample])
# [1] 415221
# [1] 109924


