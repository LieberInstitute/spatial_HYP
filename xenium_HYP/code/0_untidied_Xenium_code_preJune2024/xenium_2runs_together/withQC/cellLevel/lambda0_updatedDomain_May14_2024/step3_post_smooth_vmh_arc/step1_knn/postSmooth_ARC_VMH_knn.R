######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambdaName=0
lambda=0
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed_banksy =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"

###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")
setwd(dir_out)
dir.create("updated_May14_2024")

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
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
lambda = 0
res = 2
spe_joint = readRDS(paste0("spe_joint_2runs_banksy_lambda", lambda, "_res", res, ".rds"))
########### process the data
se = spe_joint
unique(colData(se)[,"sample_id"])
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )

######################################################################### 
## find domains
######################################################################### 
# library(mgcv)
x1 = colData(spe_joint)$x_position
x2 = colData(spe_joint)$y_position
## VMH
# normcounts = assays(spe_joint)$normcounts

### VMH
VMH_neurons = c(
  "cluster 32, VMH_Lateral_border",
  "cluster 34, VMH_4",
  "cluster 23, VMH_3",
  "cluster 36, VMH_5",
  "cluster 18, VMH_1_(Excitatory)",
  "cluster 19, VMH_2_(Excitatory)",
  "cluster 40, VMH_6_Lateral_bit"
)
# clust_M0_lam0_k50_res0.6
if_VMH = array(0, dim = ncol(spe_joint))
if_VMH[colData(spe_joint)$clust_M0_lam0_k50_res2 %in% c("34", "23", "36", "18", "19") ] = 1
length(if_VMH)



###########  smooth VMH and ARC
# fastknn:
# https://github.com/davpinto/fastknn/blob/master/README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png


library("fastknn")
library("caTools")

df_ = data.frame(x1=x1,x2=x2,if_VMH=if_VMH)
sample = sample_names[1]




#### VMH
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  library(FNN)
  knn_ = get.knn(data=data.matrix(df_[which_sample, c("x1","x2")]), 
                 k=10)
  
  knn_pred = unlist(lapply(1:length(which_sample), function(xx){
    mean(if_VMH[which_sample[knn_$nn.index[xx,]]])
  }))
  
  saveRDS(knn_pred, file = paste0("updated_May14_2024/knn_predVMH_",sample,".rds") )
  
}
knn_pred_vmh_c = array(dim=length(if_VMH))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("updated_May14_2024/knn_predVMH_",sample,".rds"))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
  knn_pred_vmh_c[which_sample] = knn_pred
}

######################################################################### 
## do the plots
######################################################################### 
plt <- function(se_, cnames, cutoff_,
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  
  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  plot_bank1 <- plotColData(se_[,grp1],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank2 <- plotColData(se_[,grp2],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  
  plot_bank_byClu1 <- plotColData(se_[,grp1],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu2 <- plotColData(se_[,grp2],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)
  
  print(plot_bank_byClu1)
  print(plot_bank_byClu2)
  
}




se = spe_joint
colDa
pdf("updated_May14_2024/multisample_banksy_2runs_13samples_domain_knn_pred_vmh_c_k10.pdf", width = 60, height = 100/2)
for(cutoff in 0.2){
  print(cutoff)
  colData(se)$pred_vmh = as.character(knn_pred_vmh_c>cutoff)

  plt(se[,], cnames = "pred_vmh",cutoff_=cutoff)
}
dev.off()



