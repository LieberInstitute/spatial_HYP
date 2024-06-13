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
dir_out = "/dcs04/hansen/data/ywang/xenium_newData/out/"
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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambda,".rds"))

se = spe_joint

######################################################################### 
## do the plots
######################################################################### 
plt <- function(){
  
  cnames <- colnames(colData(se))
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
  colData(se) <- cbind(colData(se), spatialCoords(se))
  
  plot_bank <- plotColData(se,
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu <- plotColData(se,
                                 x = "x_location", y = "y_location",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = 7)+   
    scale_shape_manual(values=seq(0,15))
  
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
  rdnames <- reducedDimNames(se)
  
  umap_ <- plotReducedDim(se,
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank)
  print(plot_bank_byClu)
  
  print(umap_)
}
# sampleID_unique  = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
pdf("multisample_banksy_lambda_1.pdf", width = 40, height = 60)
plt()
dev.off()



pdf("multisample_banksy_lambda_1_largerPage.pdf", width = 80, height = 60)
plt()
dev.off()


####
colData(se)$y_position_  = - colData(se)$y_position
plt_ <- function(){
  
  cnames <- colnames(colData(se))
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
  colData(se) <- cbind(colData(se), spatialCoords(se))
  
  plot_bank <- plotColData(se,
                           x = "x_position", y = "y_position_",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu <- plotColData(se,
                                 x = "x_position", y = "y_position_",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = 7)+   
    scale_shape_manual(values=seq(0,15))
  
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
  rdnames <- reducedDimNames(se)
  
  umap_ <- plotReducedDim(se,
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank)
  print(plot_bank_byClu)
  
  print(umap_)
}

# sampleID_unique  = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
pdf("multisample_banksy_lambda_1_.pdf", width = 40, height = 60)
plt()
dev.off()



pdf("multisample_banksy_lambda_1_largerPage_.pdf", width = 80, height = 60)
plt()
dev.off()