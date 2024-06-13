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

dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)

############ load cell type clustering result, cell level
lambda <- 0
lambdaName = "0"
se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", 2, ".rds"))


####### load manual annotation result
cellTypeAnnotations = readRDS("cellTypeAnnotations.rds")

########### 
sample_names= c(c("br6197","br5993a","br5993b","br1735a","br1735b","br6588"), c("br1225a","br1225b","br8741c",  "br8741d","br5459a","br5459b", "br8667c"))
sex_perSampleID = c(c("Male", "Female", "Female", "Male", "Male", "Female"), c("Male","Male","Female","Female","Male","Male","Female"))
names(sex_perSampleID) = sample_names

####### 
# which_f =  which(sex_perSampleID[colData(se)[,"sample_id"]]=="Female")
# which_m = which(sex_perSampleID[colData(se)[,"sample_id"]]=="Male")
sex_per_cell = sex_perSampleID[colData(se)[,"sample_id"]]
cellClu_perCell = colData(se)[,"clust_M0_lam0_k50_res2"]
cellType_perCell = paste0("cluster ", cellClu_perCell,", ", cellTypeAnnotations[as.character(cellClu_perCell)])
tb_sex_cellType = table( cellType_perCell, sex_per_cell)

write.csv(tb_sex_cellType, file = "tb_sex_cellType.csv", quote = F)

####### 
tb_sampleID_cellType = table( cellType_perCell, colData(se)[,"sample_id"])
colnames(tb_sampleID_cellType) = paste0("(", sex_perSampleID[colnames(tb_sampleID_cellType)],")", colnames(tb_sampleID_cellType))
tb_sampleID_cellType = tb_sampleID_cellType[,order(colnames(tb_sampleID_cellType))]
# tb_sampleID_cellType = tb_sampleID_cellType[,order(sex_perSampleID[colnames(tb_sampleID_cellType)])]


write.csv(tb_sampleID_cellType, file = "tb_sampleID_cellType.csv", quote = F)


