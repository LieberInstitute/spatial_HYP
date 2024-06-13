######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 1

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"


###########
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
se  = readRDS(paste0("spe_joint_banksy_lambda0_UMAP.rds"))
# spe_joint  = readRDS(paste0(dir_data_processed_banksy, "se_sample",sampleID_tmp,
                     # "_preprosssed_lambda0_UMAP.rds"))



plt <- function(se){
  
  cnames <- colnames(colData(se))
  
  # cnames
  
  cnames <- cnames[grep("^clust", cnames)][1]
  colData(se) <- cbind(colData(se), spatialCoords(se))
  
  plot_bank <- plotColData(se,
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)
  
  plot_bank_byClu <- plotColData(se,
                                 x = "x_location", y = "y_location",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = 6)
  
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
  rdnames <- reducedDimNames(se)
  
  umap_ <- plotReducedDim(se,
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 3)
  
  print(plot_bank)
  print(plot_bank_byClu)
  
  print(umap_)
}

######################################################################### 
######################################################################### 
for(res in c(1, 1.2, 2, 4, 5)){
# for(res in c(0.5, 1, 1.2, 2, 4, 5)){
    
  list_colData = readRDS(paste0("list_colData_banksy_lambda0_res", res, ".rds"))
  
  # for(sample in sample_names){
  colData(se) = rbind((list_colData[[1]]),
                      (list_colData[[2]]),
                      (list_colData[[3]]),
                      (list_colData[[4]]),
                      (list_colData[[5]]),
                      (list_colData[[6]]))    
    
  # sampleID_unique  = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
  pdf(paste0("/dcs04/hansen/data/ywang/xenium/wrapup/multi_sample_banksy/multisample_banksy_lambda0_res",res,".pdf"), width = 40, height = 60)
  plt(se)
  dev.off()
    
    
  # pdf("/dcs04/hansen/data/ywang/xenium/wrapup/multi_sample_banksy/multisample_banksy_lambda_1_largerPage.pdf", width = 80, height = 60)
  # plt()
  # dev.off()
    
  # }
}

# saveRDS(list_colData, file = paste0("list_colData_banksy_lambda1_res", res, ".rds"))
# >   saveRDS(se, file = paste0(dir_data_processed_banksy, 
                              # +                             "se_sample",sampleID_tmp,
                              # +                             "_preprosssed_lambda0_UMAP.rds"))

######################################################################### 
## do the plots
######################################################################### 



