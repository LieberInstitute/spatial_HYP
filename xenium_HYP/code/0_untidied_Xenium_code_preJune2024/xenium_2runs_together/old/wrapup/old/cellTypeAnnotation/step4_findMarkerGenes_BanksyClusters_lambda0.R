######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

lambda=0
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

###########
sample_names_old = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))
se = readRDS("spe_joint_2runs_banksy_lambda0_nucleus.rds")


######################################################################### 
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se[,colData(se)$sizeFactor>0])
se <- scater::logNormCounts(se)
se_seurat  <- as.Seurat(se) 


######################################################################### 
######################################################################### 

cnames <- colnames(colData(se))

cnames <- cnames[grep("^clust", cnames)]
print(cnames)
# clust_M0_lam0_k50_res0.6
# for(res in c(0.5, 1, 1.2, 2, 4, 5)){
  # for(res in c(0.5, 1, 1.2, 2, 4, 5)){

  
  Idents(se_seurat) = colData(se)[,cnames]
  se_seurat_markers = FindAllMarkers(se_seurat)
  
  saveRDS(se_seurat, file = paste0(
    "markerGenes_multiSampleBanksy_13samples_lambda",
    lambdaName, "_nucleus_res",0.6,".rds"))
  
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
                               "top10_markerGenes_multiSampleBanksy_13samples_lambda",
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
# }





