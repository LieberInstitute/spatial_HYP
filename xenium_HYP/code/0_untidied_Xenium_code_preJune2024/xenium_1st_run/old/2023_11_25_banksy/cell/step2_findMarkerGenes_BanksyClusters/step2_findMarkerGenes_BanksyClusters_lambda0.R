# ml conda_R/4.3
# R


###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"

dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# dir.create(dir_out)
setwd(dir_out)

###########
library(Seurat)
# library(Banksy)
library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)

library(scater)
library(cowplot)
library(ggplot2)




######################################################################### 
sampleID_unique = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
# > sampleID_unique
# [1] br6197  br5993a br5993b br1735a br1735b br6588 

ksample=1
for(ksample in 1){
  sampleID_tmp = sampleID_unique[ksample]
 
  se = readRDS(paste0(dir_data_processed_banksy, "se_sample",sampleID_tmp,"_preprosssed.rds"))
  
  
  se_seurat  <- as.Seurat(se) 
  
  Idents(se_seurat) = colData(se)[,"clust_M1_lam0_k50_res1.2"]
  se_seurat = FindAllMarkers(se_seurat)
  saveRDS(se_seurat, file = paste0(dir_data_processed_banksy, "se_seurat_sample",sampleID_tmp,"_markerGenesBanksy_lambda0.rds"))
}

# Warning messages:
# 1: Cannot add objects with duplicate keys (offending key: PC_), setting key to 'pca_m1_lam0.99_' 
# 2: Cannot add objects with duplicate keys (offending key: PC_), setting key to 'pca_m1_lam1_' 
# se_seurat = readRDS(paste0(dir_data_processed_banksy, "se_seurat_sample",sampleID_tmp,"_markerGenesBanksy.rds"))

######################################################################### 
### plot the results - heatmap for marker genes
######################################################################### 
library(dplyr)



pdf("markerGenes_top10_banksy_lambda_0.pdf", width = 20, height = 20)
for(ksample in 1:6){
  print(ksample)

  sampleID_tmp = sampleID_unique[ksample]
  se = readRDS(paste0(dir_data_processed_banksy, "se_sample",sampleID_tmp,"_preprosssed.rds"))
  
  ######## format se into seurat format
  obj = CreateSeuratObject(assays(se)$counts, project = paste0(sampleID_tmp,", marker genes"))
  # obj@data = 
  obj = NormalizeData(obj)
  # obj@data 
  obj = ScaleData(obj)
  Idents(obj) = colData(se)[,"clust_M1_lam0_k50_res1.2"]
    
  ######## 
  se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, "se_seurat_sample",sampleID_tmp,"_markerGenesBanksy_lambda0.rds"))
  
  se_seurat_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  # se_seurat  <- as.Seurat(se) 
  # VariableFeatures(se_seurat) = rownames(se)
  # se_seurat  = ScaleData(se_seurat)
  set.seed(111)#plot a subset of cells to reduce running time
  rd_cells = sample(1:ncol(obj), 10000)
  print(DoHeatmap(obj[top10$gene, rd_cells] , features = top10$gene) + NoLegend())

}
dev.off()

# se_seurat@logcounts

# se_seurat@data = assays(se)$logcounts

assays(se)

