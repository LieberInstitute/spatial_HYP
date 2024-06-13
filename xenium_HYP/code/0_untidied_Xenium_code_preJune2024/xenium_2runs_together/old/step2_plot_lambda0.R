######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 0

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
sample_names_old = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_2runs_banksy_lambda",lambda,".rds"))
# spe_joint_2runs_banksy_lambda0.rds
se = spe_joint
unique(colData(se)[,"sample_id"])
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )

######################################################################### 
## do the plots
######################################################################### 
plt <- function(se_, 
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  cnames <- colnames(colData(se_))
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
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
  
  umap_1 <- plotReducedDim(se_[,grp1],
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 2)+   
    scale_shape_manual(values=seq(0,15))
  
  umap_2<- plotReducedDim(se_[,grp2],
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 2)+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)

  print(plot_bank_byClu1)
  print(plot_bank_byClu2)
  
  print(umap_1)
  print(umap_2)
  
}
# sampleID_unique  = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
pdf("multisample_banksy_2runs_lambda_0.pdf", width = 60, height = 100/2)
plt(se[,])
# plt(se[,])
dev.off()



pdf("multisample_banksy_2runs_lambda_0_largerPage.pdf", width = 120, height = 60/2)
plt(se)
dev.off()

# 
# ####
# colData(se)$y_position_  = - colData(se)$y_position
# plt_ <- function(){
#   
#   cnames <- colnames(colData(se))
#   
#   cnames <- cnames[grep("^clust", cnames)]
#   print(cnames)
#   
#   colData(se) <- cbind(colData(se), spatialCoords(se))
#   
#   plot_bank <- plotColData(se,
#                            x = "x_position", y = "y_position_",
#                            point_size = 0.6, colour_by = cnames, shape = "sample_id"
#   ) + coord_equal() + 
#     ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
#     theme(plot.title = element_text(size=33)) + 
#     facet_wrap(~ shape, ncol = 3)+   
#     scale_shape_manual(values=seq(0,15))
#   
#   plot_bank_byClu <- plotColData(se,
#                                  x = "x_position", y = "y_position_",
#                                  point_size = 0.6, colour_by = cnames, shape = "sample_id"
#   ) + coord_equal()+ 
#     ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
#     theme(plot.title = element_text(size=33)) + 
#     facet_wrap(~ colour_by + shape, ncol = 3)+   
#     scale_shape_manual(values=seq(0,15))
#   
#   
#   # Visualize UMAPs of the non-spatial and BANKSY embedding:
#   rdnames <- reducedDimNames(se)
#   
#   umap_ <- plotReducedDim(se,
#                           dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
#                           colour_by = cnames, shape = "sample_id"
#   ) + 
#     ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
#     theme(plot.title = element_text(size=33)) + 
#     facet_wrap(~shape, ncol = 3)+   
#     scale_shape_manual(values=seq(0,15))
#   
#   print(plot_bank)
#   print(plot_bank_byClu)
#   
#   print(umap_)
# }
# 
# # sampleID_unique  = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
# pdf("multisample_banksy_2runs_lambda_0_.pdf", width = 60, height = 100)
# plt()
# dev.off()
# 
# 
# 
# pdf("multisample_banksy_2runs_lambda_0_largerPage_.pdf", width = 120, height = 60)
# plt()
# dev.off()
# 
# 
