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
########### load ARC annotation result
knn_pred_arc_c_ = readRDS("knn_ARC_twoRounds.rds")
########################################################################
########### load cell clusters within ARC
cellClu_kp = readRDS("cellClu_kp_ARC.rds")

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
which_cells_kp =  which( knn_pred_arc_c_ =="ARC")

##
cellClu_domain_kp = knn_pred_arc_c_[which_cells_kp]
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
  # which_f= sex_perCell_tmp=="Female"
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
saveRDS(list_DEG_perClu, file = "list_DEG_perClu.rds")
saveRDS(list_pval_allGenes_perClu, file = "list_pval_allGenes_perClu.rds")
saveRDS(list_fc_allGenes_perClu, file = "list_fc_allGenes_perClu.rds")
list_DEG_perClu = readRDS("list_DEG_perClu.rds")

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
saveRDS(pval_, file = "pval_allCellsWithinARC.rds")
saveRDS(fc, file = "fc_allCellsWithinARC.rds")
saveRDS(genes_deg, file = "genes_deg_allCellsWithinARC.rds")

list_DEG_perClu_cellTypeSpecific = list()
for(cellType_tmp in names(list_DEG_perClu)){
  genes_deg_tmp = list_DEG_perClu[[cellType_tmp]]
  list_DEG_perClu_cellTypeSpecific[[cellType_tmp]] = genes_deg_tmp[!genes_deg_tmp%in%genes_deg]
}

saveRDS(list_DEG_perClu_cellTypeSpecific, file = "list_DEG_perClu_cellTypeSpecific_WithinARC.rds")

# > list_DEG_perClu_cellTypeSpecific
# $`12`
# [1] "CRH"      "RELN"     "TRH"      "VGF"      "C1QL3"    "AGRP"    
# [7] "CYP26A1"  "NPY1R"    "TRPC6"    "RGS16"    "SSTR1"    "IGFBP5"  
# [13] "TRPC5"    "TMEM132C" "APOE"     "BTBD11"  
# 
# $`17`
# [1] "TAC1"   "THBS1"  "SULF1"  "PDGFRA" "NR4A2"  "IFITM3" "RGS16"  "ALKAL2"
# [9] "FGFR3"  "DNER"   "HRH1"   "EFHD1"  "SOX11"  "SNTB2" 
# 
# $`25`
# [1] "GHRH"   "VGF"    "ECEL1"  "PDGFRA" "CHODL"  "PCSK1"  "NPY1R"  "RGS16" 
# [9] "HTR2C"  "CDH12"  "APOE"   "VWC2L"  "PLD5"   "SSTR2"  "TSHZ2"  "BTBD11"
# [17] "STAT3" 
# 
# $`2`
# [1] "PDGFRA"  "CYP26A1" "IGFBP5"  "SOX10"   "VCAN"    "ERBB3"   "APOE"   
# [8] "HES1"    "EGFR"    "SHANK3"  "STAT3"  
# 
# $`1`
# character(0)
# 
# $`5`
# [1] "NR4A2"    "CAV1"     "KLF2"     "IFITM3"   "PECAM1"   "IGFBP5"  
# [7] "IGFBP4"   "TMEM132C" "KLF4"     "HES1"     "CTSH"     "MT2A"    
# [13] "KRT8"     "EFHD1"   
# 
# $`33`
# [1] "VIP"      "IGFBP3"   "RSPO2"    "THBS1"    "VGF"      "C1QL3"   
# [7] "ECEL1"    "NR4A2"    "KIT"      "CHGB"     "PCSK1"    "CYP26A1" 
# [13] "COL25A1"  "NPY1R"    "GLRA3"    "CDH13"    "IFITM3"   "TRPC6"   
# [19] "SSTR1"    "SNCA"     "HS3ST4"   "IGFBP5"   "LYPD6"    "ALK"     
# [25] "PLCH1"    "LYPD6B"   "MC3R"     "ALKAL2"   "PCDH19"   "PROKR1"  
# [31] "CNTN2"    "RASGRP1"  "SHISAL2B" "NNAT"     "CYP26B1"  "TRPC5"   
# [37] "TMEM132C" "APOE"     "AR"       "ENOX2"    "MYO5B"    "PLCE1"   
# [43] "APP"      "CDH6"     "NPNT"     "CNTNAP3B" "SLC24A3"  "PLD5"    
# [49] "TENM1"    "BRINP3"   "TRIL"     "NRP1"     "MAN1A1"   "KCNH5"   
# [55] "ANGPT1"   "HRH3"     "BTBD11"   "KRT18"    "GFOD1"    "CCNA1"   
# [61] "CDH4"     "ADRA1A"   "SOX11"    "TGFB1"    "STAT3"    "SNCG"    
# [67] "CELF2"    "ADARB1"  









