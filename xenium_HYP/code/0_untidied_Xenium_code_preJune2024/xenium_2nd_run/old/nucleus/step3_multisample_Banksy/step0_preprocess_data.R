# ml conda
# conda activate singleCellTK
# R

# following tyhis tutorial for multi-sample Banksy:
# https://rdrr.io/github/prabhakarlab/Banksy/f/vignettes/multi-sample.Rmd
###########
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"

dir_data_processed = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
dir_out = "/dcs04/hansen/data/ywang/xenium_newData/out/"
setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)

######################################################################### 
spe = readRDS(paste0(dir_data_processed, 
                       "xens_nucleus.rds")) # spe obj, not too much QC filteration
sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))

######################################################################### 
#################### step1 create spe_list
######################################################################### 
imgData(spe) <- NULL
assay(spe, "logcounts") <- NULL
reducedDims(spe) <- NULL
rowData(spe) <- NULL
# colData(spe) <- DataFrame(
#   
#   sample_id = spe$sample_id,
#   clust_annotation = factor(
#     addNA(spe$layer_guess_reordered_short),
#     exclude = NULL, labels = seq(8)
#   ),
#   in_tissue = spe$in_tissue,
#   row.names = colnames(spe)
# )
# gc()

# Next, subset spe to samples from the last subject (samples 151673, 151674, 151675, 151676). This stores each sample in its own SpatialExperiment object, all placed in a list.

# sample_names <- as.character(151673:151676)
spe_list <- lapply(sample_names, function(x) spe[, spe$sample_id == x])
# rm(spe)
# gc()

######################################################################### 
# Data preprocessing
######################################################################### 
# Using the singleCellTK package, we perform basic normalisation of the data:
######################################################################### 
library(singleCellTK)
aname <- "normcounts"
spe_list <- lapply(spe_list, function(x) {
  singleCellTK::runSeuratNormalizeData(
    x,
    useAssay = "counts",
    normAssayName = aname,
    normalizationMethod = "RC",
    scaleFactor = 5000,
    verbose = FALSE
  )
})


# Next, select the top 2000 highly variable features from each sample, and take their union for downstream analysis:
hvgs <- lapply(spe_list, function(x) {
  singleCellTK::getTopHVG(
    singleCellTK::runSeuratFindHVG(x, hvgNumber = 2000, verbose = FALSE),
    hvgNumber = 2000
  )
})
hvgs <- Reduce(union, hvgs)
# Subset to hvgs
spe_list <- lapply(spe_list, function(x) x[hvgs, ])


saveRDS(spe_list, file = paste0(dir_data_processed_banksy,"spe_list_nucleus.rds"))
