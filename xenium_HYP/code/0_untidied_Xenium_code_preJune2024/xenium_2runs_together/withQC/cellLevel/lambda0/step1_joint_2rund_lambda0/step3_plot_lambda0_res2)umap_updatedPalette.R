######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 0
res = 2
###########
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
# dir.create(dir_data_processed_banksy)

dir_data_processed = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
# dir_out = "/dcs04/hansen/data/ywang/xenium_newData/out/"
# setwd(dir_out)
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"
# dir.create(dir_data_processed_banksy)
# dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"

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
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
spe_joint = readRDS(paste0("spe_joint_2runs_banksy_lambda", lambda, "_res", res, ".rds"))
# spe_joint_2runs_banksy_lambda0.rds
se = spe_joint
unique(colData(se)[,"sample_id"])
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )

######################################################################### 
########### load manual cell type annotation info
# saveRDS(cellTypeAnnotations, file = "cellTypeAnnotations.rds")
cellTypeAnnotations = readRDS("cellTypeAnnotations.rds")
# saveRDS(vec_large_clus, file = "vec_large_clus.rds")
cellTypeAnnotations = readRDS("vec_large_clus.rds")

######################################################################### 
## plot all cells , color by large cell group
######################################################################### 
categorized_cell_types = readRDS("categorized_cell_types.rds")
cellTypeAnnotations = readRDS("cellTypeAnnotations.rds")
cellTypeAnnotations_ = cellTypeAnnotations
for(cellType_large_tmp in names(categorized_cell_types)){
  name_cellType_tnp = sub(".*, ","",categorized_cell_types[[cellType_large_tmp]])
  cellTypeAnnotations_[cellTypeAnnotations_%in% name_cellType_tnp] = cellType_large_tmp
}

# 
library(Polychrome)
P12 = createPalette(12,  c("#ff0000", "#00ff00", "#0000ff"))

cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
print(cnames)

names(P12) = unique(cellTypeAnnotations_)

pdf(paste0("multisample_banksy_2runs_lambda_0_res",res,"_annotatedCellClu_umapOnly_updatedColors.pdf"), width = 10, height = 10)
se_ = se[,colData(se)$sample_id == "br6197"]
colData(se_) <- cbind(colData(se_), spatialCoords(se_))
colData(se_)[,cnames] = paste0(cellTypeAnnotations_[colData(se_)[,cnames]] )
rdnames <- reducedDimNames(se_)

names(P12) = unique(colData(se_)[,cnames])

umap_1 <- plotReducedDim(se_,
                         dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                         colour_by = cnames, shape_by = "sample_id"
) + 
  ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
  theme(plot.title = element_text(size=33))# + 
  scale_colour_manual(values = P12)
umap_1
dev.off()

# , legend.text = element_text(size=66)
# umap_1_2 = (reducedDim(se_, "UMAP_M0_lam0"))
# # > dim(umap_1_2)
# # [1] 72010     2
# xlim_ = range(umap_1_2[,1])
# ylim_ = range(umap_1_2[,2])

###########################################################################
#### plot cells within each large cell group, colored by small cell group
###########################################################################
categorized_cell_types = readRDS("categorized_cell_types.rds")


cellTypeAnnotations = readRDS("cellTypeAnnotations.rds")
P41 = createPalette(41,  c("#ff0000", "#00ff00", "#0000ff"))

cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
print(cnames)
names(P41) = unique(cellTypeAnnotations)


se_ = se[,colData(se)$sample_id == "br6197"]
colData(se_) <- cbind(colData(se_), spatialCoords(se_))
colData(se_)[,cnames] = paste0("cluster ", colData(se_)[,cnames],", ", cellTypeAnnotations[colData(se_)[,cnames]] )
rdnames <- reducedDimNames(se_)

pdf(paste0("multisample_banksy_2runs_lambda_0_res",res,"_annotatedCellClu_umap_each_large_cellType.pdf"), width = 10, height = 10)
for(cellType_large_tmp in names(categorized_cell_types)){
  se_tmp = se_[,colData(se_)[,cnames]%in%categorized_cell_types[[cellType_large_tmp]]]
  umap_1 <- plotReducedDim(se_tmp,
                           dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                           colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0(cellType_large_tmp ))+ 
    theme(plot.title = element_text(size=33))# + 
  scale_colour_manual(P41)
  # + xlim(xlim_)+ ylim(ylim_)
  
  print(umap_1)
}

dev.off()








