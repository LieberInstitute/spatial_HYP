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

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
se = spe_joint

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
              "br1735a","br1735b","br6588")
sex_perSampleID = c("M", "F", "F", "M", "M", "F")
names(sex_perSampleID) = sampleID_unique
# table(sex)
# F      M 
# 215382 192328 

######################################################################### 
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)
cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

se_seurat  <- as.Seurat(se) 
se_seurat@meta.data$sex = sex_perSampleID[as.character(colData(se)$sample_id)]
sex = sex_perSampleID[as.character(colData(se)$sample_id)]

####### VHM DEG
se_vhm = se_seurat[,cellClu=="4"]
Idents(se_vhm) <- se_vhm@meta.data$sex
se_vhm = FindAllMarkers(se_vhm)
saveRDS(se_vhm, file = paste0(dir_data_processed_banksy, 
                              "se_vhm_seurat_sexDEG_multiSampleBanksy_lambda", 
                              lambdaName, 
                              ".rds"))
# se_vhm["ADARB1",]

####### ARC DEG
se_arc = se_seurat[,cellClu=="3"]
Idents(se_arc) <- se_arc@meta.data$sex
se_arc = FindAllMarkers(se_arc)

saveRDS(se_arc, file = paste0(dir_data_processed_banksy, 
                              "se_arc_seurat_sexDEG_multiSampleBanksy_lambda", 
                              lambdaName, ".rds"))


se_arc["ADARB1",]
######################################################################### 
### plot the results - heatmap for marker genes
######################################################################### 
library(dplyr)
######################################## 
###################### plot vhm
######################################## 
cellClu_used = "4"

pdf(paste0("DE_gender_vhm_multisampleBanksy_lambda_",lambdaName,".pdf"), width = 20, height = 20)
  ######## format se into seurat format
  obj = CreateSeuratObject(assays(se)$counts[, cellClu == cellClu_used], 
                           project = paste0("marker genes"))
  obj = NormalizeData(obj)
  obj = ScaleData(obj)
  Idents(obj) = sex[cellClu == cellClu_used]
    
  ######## 
  se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, 
                                     "se_vhm_seurat_sexDEG_multiSampleBanksy_lambda",
                                     lambdaName,".rds"))
  
  se_seurat_markers %>%
    group_by(cluster) %>%
    # dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  set.seed(111)#plot a subset of cells to reduce running time
  rd_cells = sample(1:ncol(obj), 10000)
  print(DoHeatmap(obj[top10$gene, rd_cells] , 
                  features = top10$gene) + NoLegend())

dev.off()

######################################## 
###################### plot arc
######################################## 
cellClu_used = "3"

pdf(paste0("DE_gender_arc_multisampleBanksy_lambda_",lambdaName,".pdf"), width = 20, height = 20)
######## format se into seurat format
obj = CreateSeuratObject(assays(se)$counts[, cellClu == cellClu_used], 
                         project = paste0("marker genes"))
obj = NormalizeData(obj)
obj = ScaleData(obj)
Idents(obj) = sex[ cellClu == cellClu_used]

######## 
se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, 
                                   "se_arc_seurat_sexDEG_multiSampleBanksy_lambda",
                                   lambdaName,".rds"))

se_seurat_markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


set.seed(111)#plot a subset of cells to reduce running time
rd_cells = sample(1:ncol(obj), 10000)
print(DoHeatmap(obj[top10$gene, rd_cells] , features = top10$gene) + NoLegend())


dev.off()

