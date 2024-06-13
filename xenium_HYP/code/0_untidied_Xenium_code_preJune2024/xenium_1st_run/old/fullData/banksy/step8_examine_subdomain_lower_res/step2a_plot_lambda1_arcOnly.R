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
# spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))

spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_arcOnly.rds"))
# spe_joint_banksy_lambda1_res2.rds
se = spe_joint
# Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

head(colData(se))

colData(se)$clust_M0_lam1_k50_res0.6_smooth = NULL
# colnames(colData(se))

colData(se)[,which(colnames(colData(se))=="clust_M0_lam1_k50_res0.6")[1]] = NULL # rm the layer annotations using full data


######################################################################### 
## do the plots
######################################################################### 
plt <- function(){
  
  cnames <- colnames(colData(se))
  
  # cnames
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
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

# sampleID_unique  = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
pdf(paste0("multisample_banksy_lambda_", lambdaName, "_arcOnly.pdf"), width = 40, height = 60)
plt()
dev.off()




pdf(paste0("multisample_banksy_lambda_", lambdaName, "_largerPage_arcOnly.pdf"), width = 80, height = 60)
plt()
dev.off()

pdf(paste0("multisample_banksy_lambda_", lambdaName, "_smallerPage_arcOnly.pdf"), width = 20, height = 30)
plt()
dev.off()


