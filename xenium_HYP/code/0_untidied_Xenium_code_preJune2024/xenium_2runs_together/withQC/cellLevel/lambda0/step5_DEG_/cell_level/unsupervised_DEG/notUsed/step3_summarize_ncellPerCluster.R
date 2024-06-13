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
aname <- "normcounts"

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("M", "F", "F", "M", "M", "F")
names(sex_perSampleID) = sampleID_unique

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
setwd(dir_out)
###################################################################### 
### load stuff
###################################################################### 
########### load cell type clustering using nucle onlu
cellType_clu_nuc = readRDS("cellType_clu_nuc.rds")
# saveRDS(cellType_clu_nuc, file = "cellType_clu_nuc.rds")
###########

###########
geneInfo_ = read.csv(paste0(dir_data_processed,"data/Xenium_HYP_Candidate_Genes.csv"))

########### load cell-type-domain level DE significancy result
# write.csv(mat_bin_overall, file="mat_mat_bin_overall_perCellType_VMH_overall.csv",quote=F)
# mat_bin_overall = read.csv("mat_mat_bin_overall_perCellType_VMH_overall.csv")
# head(mat_bin_overall)
###########

###################################################################### 
################################### load multi-sample Banksy result
###################################################################### 
############ load domain segmentation result
lambda <- 1
lambdaName="1"
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
cellClu_domain = colData(spe_joint)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
se_domain = spe_joint
cellIDs_se_domain = colData(se_domain)$cell_id


############ load cell type clustering result, nucleus only
lambda <- 0
lambdaName="0"
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_nucleus.rds"))
se = spe_joint
cellID_cellClustAnnotation=colData(se)$cell_id
# cellClu_annotation_kp = colData(se)[cellIDs_kp,paste0("clust_M0_lam", 0, "_k50_res0.6")]

########### keep cells that are in ARC domain and are kept in both process - full data and nucleo
cellIDs_kp = intersect(cellID_cellClustAnnotation,
                       cellIDs_se_domain[cellClu_domain%in%c("4","3")])
length(cellIDs_kp)
# 106570

########### 
cellClu_domain_kp = colData(se_domain)[cellIDs_kp,paste0("clust_M0_lam", 1, "_k50_res0.6")]
cellClu_annotation_kp = colData(se)[cellIDs_kp,paste0("clust_M0_lam", 0, "_k50_res0.6")]
# table(cellClu_domain_kp,cellClu_annotation_kp)
sampleID_kp = colData(se_domain)[cellIDs_kp,paste0("sample_id")]


cellClu_domain_kp_kp_trans = as.character(cellClu_domain_kp)
cellClu_domain_kp_kp_trans[(cellClu_domain_kp_kp_trans)=="3"]="ARC"
cellClu_domain_kp_kp_trans[cellClu_domain_kp_kp_trans=="4"]="VMH"

# names(sex_perSampleID) = sampleID_unique

sex_perSampleID_ = sex_perSampleID
sex_perSampleID_ [sex_perSampleID_=="M"] = "Male"
sex_perSampleID_ [sex_perSampleID_=="F"] = "Female"

sampleID_kp_gender = paste0(sampleID_kp," (",sex_perSampleID_[as.character(sampleID_kp)],")")
### ARC
which_ARC = which(cellClu_domain_kp_kp_trans=="ARC")
tb_ = table( cellClu_annotation_kp[which_ARC], sampleID_kp_gender[which_ARC])
rownames(tb_) = paste0("cell cluster ", rownames(tb_))
tb_ordered  = tb_[,c(1,5,2,3,4,6)]
write.csv(tb_ordered, 
          file="ncell_each_Domain_cellType_ARC.csv",quote=F)


### VMH

which_VMH = which(cellClu_domain_kp_kp_trans=="VMH")
tb_ = table( cellClu_annotation_kp[which_VMH], sampleID_kp_gender[which_VMH])
rownames(tb_) = paste0("cell cluster ", rownames(tb_))
tb_ordered  = tb_[,c(1,5,2,3,4,6)]
write.csv(tb_ordered, 
          file="ncell_each_Domain_cellType_VMH.csv",quote=F)



tb_ = table( cellClu_annotation_kp, cellClu_domain_kp_kp_trans)
rownames(tb_) = paste0("cell cluster ", rownames(tb_))


write.csv(tb_, 
          file="ncell_each_Domain_cellType.csv",quote=F)




