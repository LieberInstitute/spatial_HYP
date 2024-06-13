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

########### 
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

########### load sc-level result of DEG
list_DEG_perClu_sc = readRDS("list_DEG_perClu.rds")
list_pval_allGenes_perClu_sc = readRDS("list_pval_allGenes_perClu.rds")
list_fc_allGenes_perClu_sc = readRDS("list_fc_allGenes_perClu.rds")

########### load pseudobulk-level result of DEG
list_DEG_perClu_peseudobulk = readRDS("list_DEG_perClu_pseudobulkLevel_ARC.rds")
list_pval_allGenes_perClu_peseudobulk = readRDS("list_pval_allGenes_perClu_pseudobulkLevel_ARC.rds")
list_fc_allGenes_perClu_peseudobulk = readRDS("list_fc_allGenes_perClu_pseudobulkLevel_ARC.rds")

########### 
cellTypes_uni = names(list_fc_allGenes_perClu)

pdf("scatter_DEG_pseudobulk_sc.pdf")
for(cellType_tmp in cellTypes_uni){
  
  ###
  cor_pval = cor(list_pval_allGenes_perClu_peseudobulk[[cellType_tmp]],
                 list_pval_allGenes_perClu_sc[[cellType_tmp]], method = "spearman")
  print(cor_pval)
  cor_pval = round(cor_pval,2)
  
  ###
  cor_fc = cor(list_fc_allGenes_perClu_peseudobulk[[cellType_tmp]],
               list_fc_allGenes_perClu_sc[[cellType_tmp]], method = "spearman")
  print(cor_fc)
  cor_fc = round(cor_fc,2)
  
  ###
  plot(- log(list_pval_allGenes_perClu_peseudobulk[[cellType_tmp]]),
       - log(list_pval_allGenes_perClu_sc[[cellType_tmp]]),
       xlab = "pseudobulk", ylab = "sc",
       main = paste0("- log pval, cell cluster ", cellType_tmp,", cor = ",cor_pval))
  abline(0,1,col="blue")
  ###
  plot(list_fc_allGenes_perClu_peseudobulk[[cellType_tmp]],
       list_fc_allGenes_perClu_sc[[cellType_tmp]],
       xlab = "pseudobulk", ylab = "sc",
       main = paste0(" log fold change , cell cluster", cellType_tmp,", cor = ", cor_fc))
  abline(0,1,col="blue")
  
  
  
}
dev.off()


