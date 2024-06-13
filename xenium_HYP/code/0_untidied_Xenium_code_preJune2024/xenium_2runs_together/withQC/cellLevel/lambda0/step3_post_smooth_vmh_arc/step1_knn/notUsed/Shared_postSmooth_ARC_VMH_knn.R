
###########
lambda <- 1

###########
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"

setwd(dir_out)

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"

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
sample_names_old = readRDS(paste0("/dcs04/hansen/data/ywang/xenium/data_processed_banksy/", "sampleID_unique.rds")) # 1st run
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds")) #2nd run
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_2runs_banksy_lambda0.rds"))

se = spe_joint
unique(colData(se)[,"sample_id"])
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )

###########  ARC -- create binary variable for indicating if one cell is annotated as cell cluster 7
if_ARC = array(0, dim = ncol(spe_joint))
if_ARC[colData(spe_joint)$clust_M0_lam0_k50_res0.6 == "7"] = 1
length(if_ARC)


###########  smooth ARC domain
library(FNN)
x1 = colData(spe_joint)$x_position
x2 = colData(spe_joint)$y_position

# 1st round of KNN
df_ = data.frame(x1=x1,x2=x2,if_ARC=if_ARC)
sample = sample_names[1]
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)

  knn_ = get.knn(data=data.matrix(df_[which_sample, c("x1","x2")]), 
          k=10)

  knn_pred = unlist(lapply(1:length(which_sample), function(xx){
    mean(if_ARC[which_sample[knn_$nn.index[xx,]]])
  }))
  saveRDS(knn_pred, file = paste0("knn_predARC__",sample,"_k10.rds") )
  
}

# concatenate the knn average across samples
knn_pred_arc_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predARC__",sample,"_k10.rds"))
  knn_pred_arc_c[which_sample] = knn_pred
}

# get the smoothed domain indicator of the 1st run of knn
cutoff_1stRun = 0.2 
knn_pred_arc_1stRun = (knn_pred_arc_c > cutoff_1stRun)


# run the 2nd run of knn with k=50
df_ = data.frame(x1=x1, x2=x2, if_ARC = knn_pred_arc_c)
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  knn_ = get.knn(data=data.matrix(df_[which_sample, c("x1","x2")]), 
                 k=50)
  
  knn_pred = unlist(lapply(1:length(which_sample), function(xx){
    mean(knn_pred_arc_1stRun[which_sample[knn_$nn.index[xx,]]])
  }))
  
  saveRDS(knn_pred, file = paste0("knn_predARC__",sample,"_k50_2ndRun.rds") )
  
}
# concatenate the knn average across samples
knn_pred_arc_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predARC__",sample,"_k50_2ndRun.rds"))
  knn_pred_arc_c[which_sample] = knn_pred
}

cutoff_2ndRun = 0.3
knn_pred_arc_2ndRun = (knn_pred_arc_c > cutoff_2ndRun)
