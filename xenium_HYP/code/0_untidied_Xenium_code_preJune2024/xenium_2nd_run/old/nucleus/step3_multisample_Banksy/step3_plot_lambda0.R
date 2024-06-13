######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaname <- 0

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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambda,"_nucleus.rds"))

se = spe_joint
cnames <- colnames(colData(se))
# [1] "brnum"                    "sample_id"               
# [3] "cell_id"                  "x_position"              
# [5] "y_position"               "detecs"                  
# [7] "ngene"                    "sizeFactor"              
# [9] "slide"                    "clust_M0_lam0_k50_res0.6"
# colData(se)[,"clust_M0_lam0_k50_res0.6"]
# Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
cellType_clu_nuc = colData(se)[,"clust_M0_lam0_k50_res0.6"]
saveRDS(cellType_clu_nuc, file = "cellType_clu_nuc.rds")

######################################################################### 
## do the plots
######################################################################### 
plt <- function(){
  
  cnames <- colnames(colData(se))
  
  # cnames
  
  cnames <- cnames[grep("^clust", cnames)]
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
pdf(paste0("multisample_banksy_lambda_",lambdaname,"_nucleus.pdf"), width = 40, height = 60)
plt()
dev.off()




pdf(paste0("multisample_banksy_lambda_",lambdaname,"_largerPage_nucleus.pdf"), width = 80, height = 60)
plt()
dev.off()

