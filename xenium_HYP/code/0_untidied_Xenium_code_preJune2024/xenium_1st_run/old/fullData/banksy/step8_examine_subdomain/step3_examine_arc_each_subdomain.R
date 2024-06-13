######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 1
lambdaName = "1"

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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_arcOnly.rds"))
se = spe_joint

head(colData(se))

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

clu_unique = as.character(sort(unique(Idents(obj))))

length(clu_unique)
# 11

for(kClu in 1:length(clu_unique)){
  print(kClu)
  clu_tmp = clu_unique[kClu]
  cluster_tmp.markers = FindMarkers(obj, ident.1 = clu_tmp)
  saveRDS(cluster_tmp.markers, 
          file= paste0("cluster",clu_tmp,"_markers_arc.rds"))
  
}

for(kClu in 1:length(clu_unique)){
  clu_tmp = clu_unique[kClu]
  print(clu_tmp)
  
  cluster_tmp.markers = readRDS(paste0("cluster",clu_tmp,"_markers_arc.rds"))
  
  cluster_tmp.markers %>%
    # group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  print(top10)
}


# [1] "1"
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# CCK     0   1.841078 0.274 0.104         0
# [1] "2"
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# TAC3        0   2.754889 0.585 0.142         0
# CRH         0   1.547743 0.303 0.118         0
# RSPO2       0   1.532548 0.559 0.275         0
# IGFBP3      0   2.217724 0.696 0.307         0
# GLRA2       0   1.326983 0.945 0.693         0
# ESR1        0   1.679526 0.779 0.414         0
# SSTR1       0   1.071010 0.844 0.593         0
# RXFP1       0   1.539882 0.284 0.095         0
# CYP26A1     0   1.902531 0.330 0.135         0
# LYPD6       0   1.240141 0.455 0.257         0
# [1] "3"
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# WIF1   0.000000e+00   1.240132 0.812 0.588  0.000000e+00
# DCN    0.000000e+00   1.069944 0.562 0.328  0.000000e+00
# SOX9   0.000000e+00   1.028127 0.818 0.594  0.000000e+00
# CRYM  3.370668e-149   1.386463 0.288 0.170 1.233664e-146
# ANXA1 2.221869e-137   1.404910 0.275 0.160 8.132040e-135
# [1] "4"
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# CNDP1   0.000000e+00   2.083062 0.451 0.143  0.000000e+00
# CLDN11  0.000000e+00   1.689265 0.463 0.247  0.000000e+00
# ERMN    0.000000e+00   2.021194 0.704 0.312  0.000000e+00
# MAG     0.000000e+00   1.930770 0.298 0.081  0.000000e+00
# MOG     0.000000e+00   1.958652 0.481 0.180  0.000000e+00
# OPALIN  0.000000e+00   1.951394 0.581 0.229  0.000000e+00
# MYRF    0.000000e+00   1.887364 0.355 0.139  0.000000e+00
# UGT8    0.000000e+00   1.612567 0.461 0.239  0.000000e+00
# MOBP    0.000000e+00   1.762846 0.767 0.414  0.000000e+00
# MAL    6.485645e-250   1.260047 0.542 0.362 2.373746e-247
# [1] "5"
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# FLT1    0.000000e+00   2.467937 0.360 0.132  0.000000e+00
# PECAM1 2.117551e-280   2.053591 0.302 0.130 7.750237e-278
# PHLDB2 6.633166e-184   1.782498 0.264 0.127 2.427739e-181
# IFITM3 1.878659e-160   1.373869 0.704 0.590 6.875892e-158
# CAV1   2.319022e-154   1.496833 0.198 0.087 8.487619e-152
# DCN    5.177953e-134   1.428504 0.481 0.343 1.895131e-131
# IGFBP4 9.356649e-132   1.584242 0.252 0.136 3.424533e-129
# GPER1  3.670072e-116   1.857693 0.212 0.111 1.343246e-113
# KLF2   2.929148e-115   1.653566 0.168 0.079 1.072068e-112
# CEMIP2  2.810016e-92   1.376549 0.341 0.235  1.028466e-89
# [1] "6"
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# POMC      0.000000e+00   1.531661 0.810 0.595  0.000000e+00
# ANKRD34B 4.156919e-183   1.197676 0.470 0.297 1.521432e-180
# HRH1     7.435140e-170   1.080037 0.393 0.229 2.721261e-167
# FEZF1    1.764924e-116   1.093798 0.362 0.229 6.459623e-114
# ADCYAP1  5.412296e-112   1.136754 0.372 0.238 1.980900e-109
# [1] "7"
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# GHRH     0   2.230069 0.783 0.261         0
# GAL      0   1.935928 0.897 0.445         0
# GSX1     0   1.313052 0.397 0.157         0
# [1] "8"
# [1] p_val      avg_log2FC pct.1      pct.2      p_val_adj 
# <0 rows> (or 0-length row.names)
# [1] "9"
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# CX3CR1 1.787207e-178   1.131962 0.429 0.246 6.541179e-176
# P2RY12 2.730096e-171   1.111422 0.411 0.234 9.992151e-169
# GPR34   8.512332e-91   1.128561 0.233 0.129  3.115513e-88
# GPR183  2.272022e-48   1.171924 0.110 0.057  8.315602e-46
# [1] "10"
# [1] p_val      avg_log2FC pct.1      pct.2      p_val_adj 
# <0 rows> (or 0-length row.names)
# [1] "11"
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# ANXA1  9.491755e-121   3.010605 0.729 0.172 3.473982e-118
# SFRP2  3.880717e-108   2.984189 0.488 0.080 1.420342e-105
# NES     1.202817e-67   2.004169 0.901 0.465  4.402311e-65
# GPNMB   3.361786e-41   1.900529 0.404 0.114  1.230414e-38
# SULF1   1.745142e-37   1.714740 0.729 0.409  6.387219e-35
# SNCA    2.671788e-36   1.254098 0.911 0.715  9.778743e-34
# IGFBP5  3.269982e-34   1.254459 0.778 0.479  1.196813e-31
# DCN     2.555156e-27   1.679717 0.655 0.354  9.351872e-25
# RSPO2   2.435580e-26   1.497179 0.611 0.319  8.914223e-24
# TGFB2   7.784906e-25   1.361779 0.764 0.523  2.849276e-22
# 
