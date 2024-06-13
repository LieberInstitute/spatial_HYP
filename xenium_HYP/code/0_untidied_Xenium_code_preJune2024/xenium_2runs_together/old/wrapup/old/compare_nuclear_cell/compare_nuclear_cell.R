######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

lambda=0
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# setwd(dir_out)

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
aname <- "normcounts"

###########
# sample_names_old = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
# sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))

###########
###########
# dir_data_processed_downstream_new_nucleud = "/dcs04/hansen/data/ywang/xenium/data_processed_newRun_nuclear/"
# dir_data_processed_downstream_old_nucleud = "/dcs04/hansen/data/ywang/data_processed_newRun_nuclear/"
# dir_data_processed_downstream_new = "/dcs04/hansen/data/ywang/xenium/data_processed_newRun/"
dir_data_processed_downstream_old = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_downstream_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"

# dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
# dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
# dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
# dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"

###########
# dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

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
spe_joint_nu_old = readRDS(paste0(dir_data_processed_downstream_old, "xens_nucleus.rds"))
spe_joint_cell_old = readRDS(paste0(dir_data_processed_downstream_old, "xens.rds"))

spe_joint_nu_new = readRDS(paste0(dir_data_processed_downstream_new, "xens_nucleus.rds"))
spe_joint_cell_new = readRDS(paste0(dir_data_processed_downstream_new, "xens.rds"))
# 
###########
# colData(spe_joint_nu)
# colData(spe_joint_cell)

counts_cell_old =assays(spe_joint_cell_old)$counts
counts_cell_new =assays(spe_joint_cell_new)$counts



# dim(counts_cell)
# [1]    366 911435
#

counts_nu_old =assays(spe_joint_nu_old)$counts
counts_nu_new =assays(spe_joint_nu_new)$counts

length(intersect(colnames(counts_cell_old), colnames(counts_nu_old)))
length(intersect(colnames(counts_cell_new), colnames(counts_nu_new)))
length(intersect(colnames(counts_nu_old), colnames(counts_cell_new)))
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
# [1] 395480373
# > sum(csum_nu[cells_shared])
# [1] 112433369

sum(csum_cell[cells_shared_subsample])
sum(csum_nu[cells_shared_subsample])
# [1] 441891
# [1] 125290


