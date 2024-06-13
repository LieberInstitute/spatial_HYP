# ml conda
# conda activate singleCellTK
# R

# following tyhis tutorial for multi-sample Banksy:
# https://rdrr.io/github/prabhakarlab/Banksy/f/vignettes/multi-sample.Rmd
###########

###########
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20231019__161052__101923_KMon101923/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out_QC/"
dir_data_processed_downstream = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
setwd(dir_out)
###########
###########

library(Seurat)
# library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)



######################################################################### 
# sfe_filterd_cells = 
######################################################################### 

spe = readRDS(paste0(dir_data_processed_downstream, 
                       "xens.rds")) # spe obj, not too much QC filteration
# sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
sample_names <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

  # xennames <- c("br1225a","br1225b","br8741c",
  #               "br8741d","br5459a","br5459b",
  #               "br8667c")
# sex = c("M", "M","F", 
        # "F", "M", "M",
        # "F")
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
dir_allSamples = dir_data
list_sample = list.files(dir_allSamples)

spe_list_kp = spe_list
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  
  # sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_keep_oneSample_",sample_name,".rds"))
  cell_id_kp = readRDS(paste0(dir_data_processed_downstream, "cell_id_kp_oneSample_",sample_name,".rds"))
  # cell_id_kp = colData(sfe)$cell_id
  # saveRDS(sfe, file = paste0(dir_data_processed_downstream, "cell_id_kp_oneSample_",sample_name,".rds"))
  spe_list_kp[[ksample]] = spe_list[[ksample]][,paste0(sample_names[ksample],".",cell_id_kp)]
}
# dim(spe_list_kp[[1]])
# > dim(spe_list_kp[[1]])
# [1]   541 72010
# > 
#   > dim(spe_list[[1]])
# [1]   541 74033

######################################################################### 
# Data preprocessing
######################################################################### 
# Using the singleCellTK package, we perform basic normalisation of the data:
######################################################################### 
library(singleCellTK)
aname <- "normcounts"
spe_list_kp <- lapply(spe_list_kp, function(x) {
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
hvgs <- lapply(spe_list_kp, function(x) {
  singleCellTK::getTopHVG(
    singleCellTK::runSeuratFindHVG(x, hvgNumber = 2000, verbose = FALSE),
    hvgNumber = 2000
  )
})
hvgs <- Reduce(union, hvgs)
# Subset to hvgs
spe_list_kp <- lapply(spe_list_kp, function(x) x[hvgs, ])


saveRDS(spe_list_kp, file = paste0(dir_data_processed_downstream,"spe_list_kp.rds"))
