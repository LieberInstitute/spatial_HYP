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
# dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
dir_data_processed =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"

setwd(dir_out)

########################################################################
########### load VMH annotation result
knn_pred_VMH_c_ = readRDS("knn_VMH_twoRounds.rds")

###################################################################### 
################################### load multi-sample Banksy result
###################################################################### 
############ load cell type clustering result, nucleus only
se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", 2, ".rds"))

cellClu_annotation = colData(se)[,"clust_M0_lam0_k50_res2"]

########################################################################
################################# format data
######################################################################### 
se <- scater::computeLibraryFactors(se)
##
which_cells_kp =  which(knn_pred_VMH_c_ =="VMH")

##
cellClu_domain_kp = knn_pred_VMH_c_[which_cells_kp]
cellClu_annotation_kp =cellClu_annotation[which_cells_kp]
sort(table(cellClu_annotation_kp))


tb_ncell_perClu = table(cellClu_annotation_kp)

cellClu_kp = names(sort(-tb_ncell_perClu))[1:20]
# > cellClu_kp
# [1] "12" "17" "25" "2"  "1"  "5"  "33"
saveRDS(cellClu_kp, file = "cellClu_kp_VMH.rds")




