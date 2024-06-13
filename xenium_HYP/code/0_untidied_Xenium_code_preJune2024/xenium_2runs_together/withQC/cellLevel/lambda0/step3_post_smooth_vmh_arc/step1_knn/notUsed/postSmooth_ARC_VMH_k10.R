######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 1

###########
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)

dir_data_processed = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
# dir_out = "/dcs04/hansen/data/ywang/xenium_newData/out/"
# setwd(dir_out)
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
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
# sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
###########
sample_names_old = readRDS(paste0("/dcs04/hansen/data/ywang/xenium/data_processed_banksy/", "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_2runs_banksy_lambda0.rds"))

# spe_joint_2runs_banksy_lambda0.rds
se = spe_joint
unique(colData(se)[,"sample_id"])
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )

######################################################################### 
## find domains
######################################################################### 
library(mgcv)
x1 = colData(spe_joint)$x_position
x2 = colData(spe_joint)$y_position
## VMH
# normcounts = assays(spe_joint)$normcounts

### VMH
# clust_M0_lam0_k50_res0.6
if_VMH = array(0, dim = ncol(spe_joint))
if_VMH[colData(spe_joint)$clust_M0_lam0_k50_res0.6 == "7"] = 1
length(if_VMH)


### ARC
# clust_M0_lam0_k50_res0.6
if_ARC = array(0, dim = ncol(spe_joint))
if_ARC[colData(spe_joint)$clust_M0_lam0_k50_res0.6 == "10"] = 1
length(if_ARC)
# length(exp_ARC)
# 911435


###########  smooth VMH and ARC
df_ = data.frame(x1=x1,x2=x2,if_VMH=if_VMH,if_ARC=if_ARC)
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  # gam: https://noamross.github.io/gams-in-r-course/chapter3
  fit_gam_arc = gam(if_ARC ~ s(x1, x2),
                    data = df_[which_sample,], # ~5min to fit
                    method = "REML")
  saveRDS(fit_gam_arc, file = paste0("fit_gam_if_ARC_sample",sample,".rds") )
  pred_arc = predict.gam(fit_gam_arc) # ~1min to predict
  saveRDS(pred_arc, file = paste0("pred_gam_if_ARC_sample",sample,".rds") )
  
  
  fit_gam_vmh = gam(if_VMH ~ s(x1, x2),
                    data = df_[which_sample,], # ~5min to fit
                    method = "REML")
  saveRDS(fit_gam_vmh, file = paste0("fit_gam_if_VMH_sample",sample,".rds") )
  pred_vmh = predict.gam(fit_gam_arc) # ~1min to predict
  saveRDS(pred_vmh, file = paste0("pred_gam_if_VMH_sample",sample,".rds") )
  
}


############################################ 
############ process the smoothed result
############################################ 
pred_arc_c =pred_vmh_c= array(dim=ncol(spe_joint))
for(sample in sample_names){
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  pred_arc  =readRDS( paste0("pred_gam_arc_sample",sample,".rds") )
  pred_arc_c[which_sample] = pred_arc
  pred_vmh  =readRDS( paste0("pred_gam_vmh_sample",sample,".rds") )
  pred_vmh_c[which_sample] = pred_vmh
  
}

set.seed(111)
if_ARC_smoothed = kmeans(pred_arc_c, 2)$cluster
colData(spe_joint)$if_ARC_smoothed = if_ARC_smoothed
set.seed(111)
if_VMH_smoothed = kmeans(pred_vmh_c, 2)$cluster
colData(spe_joint)$if_VMH_smoothed = if_VMH_smoothed
saveRDS(if_ARC_smoothed, file = paste0("domain_if_ARC_smoothed.rds") )
saveRDS(if_VMH_smoothed, file = paste0("domain_if_VMH_smoothed.rds") )




######################################################################### 
## do the plots
######################################################################### 
plt <- function(se_, cnames,
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  

  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  plot_bank1 <- plotColData(se_[,grp1],
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank2 <- plotColData(se_[,grp2],
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  
  plot_bank_byClu1 <- plotColData(se_[,grp1],
                                 x = "x_location", y = "y_location",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu2 <- plotColData(se_[,grp2],
                                 x = "x_location", y = "y_location",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
    scale_shape_manual(values=seq(0,15))
  
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
  rdnames <- reducedDimNames(se_)
  
  # umap_1 <- plotReducedDim(se_[,grp1],
  #                         dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
  #                         colour_by = cnames, shape_by = "sample_id"
  # ) + 
  #   ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
  #   theme(plot.title = element_text(size=33)) + 
  #   facet_wrap(~shape_by, ncol = 2)+   
  #   scale_shape_manual(values=seq(0,15))
  # 
  # umap_2<- plotReducedDim(se_[,grp2],
  #                         dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
  #                         colour_by = cnames, shape_by = "sample_id"
  # ) + 
  #   ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
  #   theme(plot.title = element_text(size=33)) + 
  #   facet_wrap(~shape_by, ncol = 2)+   
  #   scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)

  print(plot_bank_byClu1)
  print(plot_bank_byClu2)
  
  # print(umap_1)
  # print(umap_2)
  
}



se = spe_joint
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])


pdf("multisample_banksy_2runs_13samples_domain_if_ARC_smoothed.pdf", width = 60, height = 100/2)
plt(se[,], cnames = "if_ARC_smoothed")
dev.off()

pdf("multisample_banksy_2runs_13samples_domain_if_VMH_smoothed.pdf", width = 60, height = 100/2)
plt(se[,], cnames = "if_VMH_smoothed")
dev.off()




