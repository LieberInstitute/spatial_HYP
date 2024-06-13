# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaName="0"

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/wrapup/multi_sample_banksy/"
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


########### load multi-sample Banksy result
se  = readRDS(paste0("spe_joint_banksy_lambda0_nucleus_UMAP.rds"))

# se = spe_joint

################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)

se_seurat  <- as.Seurat(se) 


######################################################################### 
######################################################################### 
for(res in c(0.5, 1, 1.2, 2, 4, 5)){
  # for(res in c(0.5, 1, 1.2, 2, 4, 5)){
  print(res)
  
  list_colData = readRDS(paste0("list_colData_banksy_lambda0_nucleus_res", res, ".rds"))
  # saveRDS(list_colData, file = paste0("list_colData_banksy_lambda0_nucleus_res", res, ".rds"))
  
  # for(sample in sample_names){
  colData_ = rbind((list_colData[[1]]),
                      (list_colData[[2]]),
                      (list_colData[[3]]),
                      (list_colData[[4]]),
                      (list_colData[[5]]),
                      (list_colData[[6]]))    

  
  Idents(se_seurat) = colData_[,paste0("clust_M0_lam", lambda, "_k50_res", res)]
  se_seurat_markers = FindAllMarkers(se_seurat)
  
  saveRDS(se_seurat, file = paste0(
    "se_seurat_markerGenes_multiSampleBanksy_lambda",
    lambdaName, "_nucleus_res",res,".rds"))
  
  # saveRDS(se_seurat, file = paste0(dir_data_processed_banksy, 
  #                                  "se_seurat_markerGenes_multiSampleBanksy_lambda", 
  #                                  lambdaName, "_nucleus.rds"))
  
  ######################################################################### 
  ### plot the results - heatmap for marker genes
  ######################################################################### 
  library(dplyr)
  
  pdf(paste0("markerGenes_top10_multisampleBanksy_lambda_",
             lambdaName,
             "_nucleus_res",
             res,".pdf"), width = 20, height = 20)
  
  ######## format se into seurat format
  obj = CreateSeuratObject(assays(se)$counts, project = paste0("marker genes"))
  # obj@data = 
  obj = NormalizeData(obj)
  # obj@data 
  obj = ScaleData(obj)
  Idents(obj) = colData_[,paste0("clust_M0_lam", lambda, "_k50_res",res)]
  
  ######## 
  # se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda",lambdaName,"_nucleus.rds"))
  
  
  se_seurat_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  write.csv(top10, file=paste0(
                               "top10_markerGenes_multiSampleBanksy_lambda",
                               lambdaName,
                               "_nucleus_res", res, ".csv"))
  
  
  set.seed(111) #plot a subset of cells to reduce running time
  rd_cells = c()
  
  clus_unique = unique(Idents(se_seurat))
  for(clus in clus_unique){
    which_ = which(Idents(se_seurat) == clus)
    rd_cells = c(rd_cells, 
                 sample(which_, round(length(which_)/10+1) ))
  }
  
  print(DoHeatmap(obj[top10$gene, rd_cells] , features = top10$gene) + NoLegend())
  
  
  dev.off()
  
  
  # }
}





