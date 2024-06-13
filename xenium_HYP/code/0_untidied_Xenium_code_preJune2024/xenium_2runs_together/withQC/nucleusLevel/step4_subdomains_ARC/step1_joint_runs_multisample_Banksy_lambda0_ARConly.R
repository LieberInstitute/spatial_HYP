######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# setwd(dir_out)
###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
aname <- "normcounts"
library(scater)

###########
sample_names_old = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))


spe_joint = readRDS("spe_joint_2runs_banksy_lambda0.rds")
# saveRDS(spe_joint, file = paste0("spe_joint_2runs_banksy_lambda0.rds"))

############################################
########### load ARC domain info
############################################
knn_pred_arc_c = array(dim=length(if_ARC))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predARC_",sample,"_k50.rds"))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
  knn_pred_arc_c[which_sample] = knn_pred
}

if_arc = knn_pred_arc_c > 0.1


spe_joint_arc = spe_joint[,which(if_arc)]

############################################
###########
############################################
compute_agf <- FALSE
k_geom <- 6
# spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)

# We then merge the samples to perform joint dimensional reduction and clustering:
# spe_joint <- do.call(cbind, spe_list)
# rm(spe_list)
# gc()

print("11")

# When running multi-sample BANKSY PCA, the group argument may be provided. This specifies the grouping variable for the cells or spots across the samples. Features belonging to cells or spots corresponding to each level of the grouping variable will be z-scaled separately. In this case, sample_id in colData(spe_joint) gives the grouping based on the sample of origin.
lambda <- 0
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint_arc, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)

# Run UMAP on the BANKSY embedding:
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
saveRDS(spe_joint, file = paste0("spe_joint_2runs_banksy_lambda0_ARConly_.rds"))
print("22")

# Finally, obtain cluster labels for spots across all 4 samples.
res <- 0.6 

spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000) #4min per run
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
saveRDS(spe_joint, file = paste0("spe_joint_2runs_banksy_lambda0_ARConly_.rds"))


spe_joint  = readRDS("spe_joint_2runs_banksy_lambda0_ARConly_.rds")

for(res in c(  2, 1, 3, 4, 0.6)){
  print(res)
  spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000) #4min per run
  colData_ = colData(spe_joint)
  saveRDS(colData_, file = paste0("colData_spe_joint_2runs_banksy_lambda0_ARConly_res",res,".rds"))
}

#############
##### plot
#############
plt <- function(se_, cnames, res_,
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  
  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  plot_bank1 <- plotColData(se_[,grp1],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", res = ",res_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank2 <- plotColData(se_[,grp2],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", res = ",res_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))
  
  
  plot_bank_byClu1 <- plotColData(se_[,grp1],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", res = ",res_)) + 
    theme(plot.title = element_text(size=33)) + 
    # facet_grid( vars(as.character(cnames)) ,vars("sample_id"))+#, ncol =  length(sample_names_new))+
    facet_grid( vars(as.character(cnames)) ,vars("sample_id"))+#, ncol =  length(sample_names_new))+
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu2 <- plotColData(se_[,grp2],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", res = ",res_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)
  
  print(plot_bank_byClu1)
  print(plot_bank_byClu2)
  
}


se = spe_joint
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])

pdf("multisample_banksy_2runs_13samples_subdomain_arcOnly.pdf", width = 60, height = 100/2)
for(res in c( 0.6, 1, 2, 3, 4)){
  print(res)
# for(res in c(  0.6)){
  # saveRDS(colData_, file = paste0("colData_spe_joint_2runs_banksy_lambda0_ARConly_res",res,".rds"))
  colData_ = readRDS( paste0("colData_spe_joint_2runs_banksy_lambda0_ARConly_res",res,".rds"))
  print(res)
  colData(se) = colData_
  colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])

  plt(se[,], cnames = paste0("clust_M0_lam0_k50_res", res), res=res)
}
dev.off()



