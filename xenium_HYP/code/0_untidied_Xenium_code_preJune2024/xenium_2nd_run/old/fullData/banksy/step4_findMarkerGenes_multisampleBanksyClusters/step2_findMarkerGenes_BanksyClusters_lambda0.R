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


########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))

se = spe_joint
# Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

######################################################################### 
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)

se_seurat  <- as.Seurat(se) 
Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
se_seurat = FindAllMarkers(se_seurat)

saveRDS(se_seurat, file = paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, ".rds"))

######################################################################### 
### plot the results - heatmap for marker genes
######################################################################### 
library(dplyr)
library(ggplot2)
se_seurat = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, ".rds"))

pdf(paste0("markerGenes_top10_multisampleBanksy_lambda_",lambdaName,".pdf"), width = 20, height = 20)

  ######## format se into seurat format
  obj = CreateSeuratObject(assays(se)$counts, project = paste0("marker genes"))
  # obj@data = 
  obj = NormalizeData(obj)
  # obj@data 
  obj = ScaleData(obj)
  Idents(obj) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
    
  ######## 
  # require(devtools)
  # install_version("ggplot2", version = "3.4.4")
  se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda",lambdaName,".rds"))
  
  se_seurat_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  

  set.seed(111)#plot a subset of cells to reduce running time
  rd_cells = sample(1:ncol(obj), 10000)
  print(DoHeatmap(obj[top10$gene, rd_cells] , features = top10$gene) + NoLegend())

dev.off()

### based on the spatial location of clu 13, it should be corresponding to clu15 in the 1st run
pdf(paste0("markerGenes_top10_multisampleBanksy_lambda_",lambdaName,"_check_clu13.pdf"), width = 20, height = 20)
#plot a subset of cells to reduce running time
set.seed(111)
rd_cells = sample(1:ncol(obj), 10000)
clu13_cells = which(Idents(obj)  == "13")
# set.seed(111)
clu13_cells_sub = sample(clu13_cells, 1000)
rd_cells=c(clu13_cells_sub,rd_cells)
rd_cells = rd_cells[!duplicated(rd_cells)]
genes_interested = toupper(c("gad1", "gad2", "slc32a1", "slc17a6", "slc17a7"))
genes_interested_kp = intersect(genes_interested, rownames(obj))
clu_uni = as.character( unique(Idents(obj)))
Idents(obj) = factor(as.character(Idents(obj)), levels = c("13", sort(clu_uni[clu_uni!="13"])))
print(DoHeatmap(obj[c(genes_interested, top10$gene), rd_cells] , features = c(genes_interested, top10$gene)) + NoLegend())
dev.off()


# check gad1 and slc17a6 in clu13
SLC17A6_clu13 = (obj[["RNA"]]@data["SLC17A6", clu13_cells])
summary(SLC17A6_clu13)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   0.000   1.656   3.497   6.225 
GAD1_clu13 = (obj[["RNA"]]@data["GAD1", clu13_cells])
summary(GAD1_clu13)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   4.328   5.012   4.529   5.435   6.772 
sum(SLC17A6_clu13>0 & GAD1_clu13>0)
sum(SLC17A6_clu13>0 & GAD1_clu13==0)
sum(SLC17A6_clu13==0 & GAD1_clu13>0)
sum(SLC17A6_clu13==0 & GAD1_clu13==0)
# [1] 5676
# [1] 589
# [1] 7455
# [1] 572
length(clu13_cells)
# 14292

# toupper("slc17a6")%in%rownames(obj)
# Cell clu13 in the 2nd-run data should be corresponding to the cluster 15 (the ARC-specific cell cluster with many cells expressing both GAD1 and SLC17A6) in the 1st-run data, according to its spatial locations (I'll need to jointly analyze both runts to find the exactly matched cell types)
# within cell cluster 13 (14292 cells in total) in the 2nd run data, there are 5676 cells expressing both GAD1 and SLC17A6




