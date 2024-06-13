# ml conda_R/4.3
# R
###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(cowplot)
########### 
# lambdaName=0

aname <- "normcounts"
########### 
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names= c(sample_names_old, sample_names_new)
sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female",
                    c("Male","Male","Female","Female","Male","Male","Female"))
names(sex_perSampleID) = sample_names

###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"

###########
setwd(dir_out)

########################################################################
########### load VMH annotation result
knn_pred_VMH_c_ = readRDS("knn_VMH_twoRounds.rds")
########################################################################
########### load cell clusters within VMH
cellClu_kp = readRDS("cellClu_kp_VMH.rds")

########################################################################
###################################################################### 
################################### load multi-sample Banksy result
###################################################################### 
############ load cell type clustering result, nucleus only

se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res2.rds"))

cellClu_annotation = colData(se)[,"clust_M0_lam0_k50_res2"]

########################################################################
################################# format data
######################################################################### 
se <- scater::computeLibraryFactors(se)
##
which_cells_kp =  which( knn_pred_VMH_c_ =="VMH")

##
cellClu_domain_kp = knn_pred_VMH_c_[which_cells_kp]
cellClu_annotation_kp = as.character(cellClu_annotation[which_cells_kp])
##
se <- scater::logNormCounts(se[,which_cells_kp])

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]

######## get DEG - unsupervised
list_DEG_perClu = list()
list_pval_allGenes_perClu = list()
list_fc_allGenes_perClu = list()
cell_clu = "33"
for(cell_clu in cellClu_kp){
  
  print(cell_clu)
  
  spe_tmp = se[, cellClu_annotation_kp==cell_clu]
  sex_perCell_tmp = sex_perSampleID[colData(spe_tmp)$sample_id]
  which_f= sex_perCell_tmp=="Female"
  which_m= sex_perCell_tmp=="Male"
  
  
  logcounts_spe_tmp = assays(spe_tmp)$logcounts
  logcounts_spe_tmp = assays(spe_tmp)$logcounts
  
  pval_ = unlist(lapply(rownames(logcounts_spe_tmp),function(xx){
    t.test(logcounts_spe_tmp[xx,which_m],logcounts_spe_tmp[xx,which_f])$p.value
  }))
  
  fc = unlist(lapply(rownames(logcounts_spe_tmp),function(xx){
    mean(logcounts_spe_tmp[xx,which_m])-mean(logcounts_spe_tmp[xx,which_f])
  }))
  
  genes_deg = rownames(logcounts_spe_tmp)[pval_<.05 & abs(fc)>.2]
  
  # save DEG list for each cell cluster
  list_DEG_perClu[[as.character(cell_clu)]] = genes_deg
  
  #save pval and fc for all genes
  names(pval_) = names(fc) = rownames(logcounts_spe_tmp)
  list_pval_allGenes_perClu[[as.character(cell_clu)]] = pval_
  list_fc_allGenes_perClu[[as.character(cell_clu)]] = fc
  
  print(length(genes_deg))
}
saveRDS(list_DEG_perClu, file = "list_DEG_perClu_VMH.rds")
saveRDS(list_pval_allGenes_perClu, file = "list_pval_allGenes_perClu_VMH.rds")
saveRDS(list_fc_allGenes_perClu, file = "list_fc_allGenes_perClu_VMH.rds")
list_DEG_perClu = readRDS("list_DEG_perClu_VMH.rds")

# cp list_DEG_perClu.rds 

# cp list_pval_allGenes_perClu.rds #
# cp list_fc_allGenes_perClu.rds  
#log fold change



# [1] "12"
# [1] 39
# [1] "17"
# [1] 47
# [1] "25"
# [1] 31
# [1] "2"
# [1] 35
# [1] "1"
# [1] 14
# [1] "5"
# [1] 40
# [1] "33"
# [1] 88


######## get DEG - find genes that are only cell-level DEG


# print(cell_clu)
spe_tmp = se[, ]
sex_perCell_tmp = sex_perSampleID[colData(spe_tmp)$sample_id]
which_f= sex_perCell_tmp=="Female"
which_m= sex_perCell_tmp=="Male"

logcounts_spe_tmp = assays(spe_tmp)$logcounts
logcounts_spe_tmp = assays(spe_tmp)$logcounts

pval_ = unlist(lapply(rownames(logcounts_spe_tmp),function(xx){
  t.test(logcounts_spe_tmp[xx,which_m],logcounts_spe_tmp[xx,which_f])$p.value
}))

fc = unlist(lapply(rownames(logcounts_spe_tmp),function(xx){
  mean(logcounts_spe_tmp[xx,which_m])-mean(logcounts_spe_tmp[xx,which_f])
}))

genes_deg = rownames(logcounts_spe_tmp)[pval_<.05 & abs(fc)>.2]

names(pval_) = names(fc) = rownames(logcounts_spe_tmp)

print(length(genes_deg))
saveRDS(pval_, file = "pval_allCellsWithinVMH.rds")
saveRDS(fc, file = "fc_allCellsWithinVMH.rds")
saveRDS(genes_deg, file = "genes_deg_allCellsWithinVMH.rds")

list_DEG_perClu_cellTypeSpecific = list()
for(cellType_tmp in names(list_DEG_perClu)){
  genes_deg_tmp = list_DEG_perClu[[cellType_tmp]]
  list_DEG_perClu_cellTypeSpecific[[cellType_tmp]] = genes_deg_tmp[!genes_deg_tmp%in%genes_deg]
}

saveRDS(list_DEG_perClu_cellTypeSpecific, file = "list_DEG_perClu_cellTypeSpecific_WithinVMH.rds")
