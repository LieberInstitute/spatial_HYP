######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 1
lambdaName="1"
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)
aname <- "normcounts"

###########
sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_vmhOnly.rds"))
se = spe_joint

head(colData(se))

# colData(se)$clust_M0_lam1_k50_res0.6_smooth = NULL

colData(se)[,which(colnames(colData(se))=="clust_M0_lam1_k50_res0.6")[1]] = NULL # rm the layer annotations using full data
# 
######################################################################### 
## examine DEGs for subdomain1 of vmh
######################################################################### 
library(dplyr)

obj = CreateSeuratObject(assays(se)$counts)
# obj@data = 
obj = NormalizeData(obj)
# obj@data 
obj = ScaleData(obj)
Idents(obj) = colData(se)[,"clust_M0_lam1_k50_res0.6"]
cluster1.markers = FindMarkers(obj, ident.1 = "1")

saveRDS(cluster1.markers, 
        file= "cluster1_markers_vmh.rds")


se_seurat_markers = cluster1.markers
se_seurat_markers %>%
  # group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  # slice_head(n = 10) %>%
  ungroup() -> top10


# p_val avg_log2FC pct.1 pct.2     p_val_adj
# GABRQ  0.000000e+00   1.021548 0.668 0.489  0.000000e+00
# SST   2.111478e-282   1.203820 0.426 0.271 7.728011e-280
# GABRE 2.313727e-277   1.328868 0.305 0.166 8.468242e-275
# VGF   9.633332e-245   1.163724 0.367 0.228 3.525800e-242



# p_val avg_log2FC pct.1 pct.2     p_val_adj
# GABRQ  0.000000e+00   1.021548 0.668 0.489  0.000000e+00
# SST   2.111478e-282   1.203820 0.426 0.271 7.728011e-280
# GABRE 2.313727e-277   1.328868 0.305 0.166 8.468242e-275
# VGF   9.633332e-245   1.163724 0.367 0.228 3.525800e-242

