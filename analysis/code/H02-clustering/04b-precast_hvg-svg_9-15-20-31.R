library(data.table)
library(Biostrings)
library(Seurat)
library(scater) # addPerCellQC
library(scran)
library(PRECAST)

#### Load list of 10 and 20%ile HVGs, SVGs; load list of seurat-formatted samplewise data
#### before seurat conversion, spe colnames and rowData(colnames) were changed to spe$key to ensure unique values when rejoining
#### before seurat conversion, coldata added under new entries 'row' and 'col' corresponding to $array_row and $array_col

feats <- readRDS("hvg_svg_sets_list.RDS")
seulist <- load("seurat_list_forprecast.RData")

### build table with constants (k value from among the optimized values of 9/15/20/31 and feature set to retrieve)

f <- rep(c(1:4),4)
k <- c(rep(9,4),rep(15,4),rep(20,4),rep(31,4))
fks <- as.data.frame(cbind(f,k))
rm(f,k)


### pull params for task number

i <- as.numeric(Sys.getenv("SGE_TASK_ID"))

### And using PRECAST, with more of tony's code: https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/code/EDA/01_PRECAST.R
### as well as Erik's code for bypassing SPARK-X and passing in other genes
# providing a custom gene list suppresses SPARK-X.

preobj = CreatePRECASTObject(seulist,customGenelist = featlist[[fks[i,1]]])

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 8, maxIter = 30, verbose = TRUE)


PRECASTObj <- PRECAST(PRECASTObj, K = fks[i,2])
save(PRECASTobj,paste0("precast_hvg-svg_out/",names(feats)[fks[i,1]],"_k",fks[i,2],"_precast_clusts.RData"))

rm(list=ls())
gc(full=T)
session.info()
