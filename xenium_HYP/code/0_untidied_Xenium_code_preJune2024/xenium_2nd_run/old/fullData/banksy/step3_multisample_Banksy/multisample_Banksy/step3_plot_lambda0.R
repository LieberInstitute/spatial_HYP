######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaName="0"
###########
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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
# spe_joint_banksy_lambda1_res2.rds
se = spe_joint

head(colData(se)[colData(se)$sample_id=="br8741d",])
table(colData(se)$sample_id)

table(colData(se)[colData(se)$sample_id=="br8741d","clust_M0_lam1_k50_res0.6" ])
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 705  5408  1637   190 14208  7813  2039   709    50 13705    10  2093  2530 
# 14    15    16 
# 2477  2229    49 

# br1225a br1225b br5459a br5459b br8667c br8741c br8741d 
# 97379   72855   69815   70514   74522   62788   55852 

# colData(se)$clust_M0_lam1_k50_res0.6 = NULL

######################################################################### 
## do the plots
######################################################################### 
plt <- function(){
  
  cnames <- colnames(colData(se))
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
  colData(se) <- cbind(colData(se), spatialCoords(se))
  
  plot_bank <- plotColData(se,
                           x = "x_position", y = "y_position",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu <- plotColData(se,
                                 x = "x_position", y = "y_position",
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
pdf(paste0("multisample_banksy_lambda_", lambdaName, ".pdf"), width = 40, height = 60)
plt()
dev.off()




pdf(paste0("multisample_banksy_lambda_0_largerPage", lambdaName, ".pdf"), width = 80, height = 60)
plt()
dev.off()

pdf(paste0("multisample_banksy_lambda_0_smallerPage", lambdaName, ".pdf"), width = 20, height = 30)
plt()
dev.off()


#####
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
pdf(paste0("multisample_banksy_lambda_", lambdaName, "_.pdf"), width = 40, height = 60)
plt_()
dev.off()

pdf(paste0("multisample_banksy_lambda_0_largerPage", lambdaName, "_.pdf"), width = 80, height = 60)
plt_()
dev.off()

pdf(paste0("multisample_banksy_lambda_0_smallerPage", lambdaName, "_.pdf"), width = 20, height = 30)
plt_()
dev.off()


