# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaName="0"

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


########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_nucleus.rds"))

se = spe_joint

################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)

se_seurat  <- as.Seurat(se) 
Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
se_seurat = FindAllMarkers(se_seurat)

saveRDS(se_seurat, file = paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, "_nucleus.rds"))

se_seurat = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, "_nucleus.rds"))
######################################################################### 
### plot the results - heatmap for marker genes
######################################################################### 
library(dplyr)
genes_interest = c("gad1", "gad2", "slc32a1", "slc17a6",  "slc17a7")
genes_interest = toupper(genes_interest)

########
pdf(paste0("multisampleBanksy_lambda_",lambdaName,"_nucleus_cluster15Examine.pdf"), width = 20, height = 20)
obj = CreateSeuratObject(assays(se)$counts, project = paste0("marker genes"))
obj = NormalizeData(obj)
obj = ScaleData(obj)
Idents(obj) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
Idents(obj) = factor(Idents(obj), levels= as.character(c(15,1:14,16:19)))

genes_interest = intersect(genes_interest, rownames(obj))

se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda",lambdaName,"_nucleus.rds"))

se_seurat_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


set.seed(111)#plot a subset of cells to reduce running time
rd_cells = sample(1:ncol(obj), 20000)
rd_cells = unique(c(rd_cells,(which(Idents(obj)=="15"))))
features_plot = c(genes_interest,(top10$gene)[!(top10$gene)%in%genes_interest])
print(DoHeatmap(obj[features_plot, rd_cells] , features =features_plot) + NoLegend())

dev.off()


