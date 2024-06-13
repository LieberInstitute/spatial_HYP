# ml conda
# conda activate singleCellTK
# R

# following tyhis tutorial for multi-sample Banksy:
# https://rdrr.io/github/prabhakarlab/Banksy/f/vignettes/multi-sample.Rmd
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"

dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"

setwd(dir_out)

###########
library(Seurat)
# library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
# library(scuttle)

# library(scater)
# library(cowplot)
# library(ggplot2)


######################################################################### 

spe = readRDS(paste0(dir_data_processed, 
                     "xens2_nucleus.rds")) # spe obj, not too much QC filteration
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
rm(spe)
gc()

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

######################################################################### 
######################################################################### 
######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R
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
aname <- "normcounts"
###########
sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))

# To run BANKSY across multiple samples, we first compute the BANKSY neighborhood feature matrices for each sample separately. We use k_geom=6 corresponding to the first-order neighbors in 10x Visium assays (k_geom=18 corresponding to first and second-order neighbors may also be used).
spe_list = readRDS(paste0(dir_data_processed_banksy,"spe_list_nucleus.rds"))

###########
compute_agf <- FALSE
k_geom <- 6
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)

# We then merge the samples to perform joint dimensional reduction and clustering:
spe_joint <- do.call(cbind, spe_list)
# rm(spe_list)
# gc()


# When running multi-sample BANKSY PCA, the group argument may be provided. This specifies the grouping variable for the cells or spots across the samples. Features belonging to cells or spots corresponding to each level of the grouping variable will be z-scaled separately. In this case, sample_id in colData(spe_joint) gives the grouping based on the sample of origin.
lambda <- 1
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)

# Run UMAP on the BANKSY embedding:
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda1_nucleus.rds"))

# Finally, obtain cluster labels for spots across all 4 samples.
res <- 0.6
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda1_nucleus.rds"))

# To compare clusters visually in the next section, run connectClusters:
# spe_joint <- connectClusters(spe_joint) # skipped this line, due to error messages 
# connectClusters: Relabel cluster labels across parameter runs to maximise their similarity.
# NA --> clust_M0_lam0.2_k50_res0.6
# Error in `[.data.frame`(clust_df, , curr_many) : 
  # undefined columns selected


# Once joint clustering is performed, we split the samples into their own SpatialExperiment objects:
spe_list <- lapply(sample_names, function(x) spe_joint[, spe_joint$sample_id == x])
# rm(spe_joint)
# gc()
# 
# As an optional step, we smooth the cluster labels of each sample separately. This can be done if smooth spatial domains are expected in the biological sample or tissue in question.
spe_list <- lapply(spe_list, smoothLabels, cluster_names = cnm, k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_names)

saveRDS(spe_list, file = paste0("spe_list_banksy_lambda1_nucleus.rds"))

######################################################################### 
# Parsing BANKSY output
######################################################################### 
spe_list = readRDS("spe_list_banksy_lambda1_nucleus.rds")
use_agf <- FALSE
lambda <- 1
res <- 0.6
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)

# We can compare BANKSY clusters to pathology annotations using several cluster comparison measures such as the adjusted Rand index (ARI) or normalized mutual information (NMI) with compareClusters.
# 
# ari <- sapply(spe_list, function(x) as.numeric(tail(compareClusters(x, func = "ARI")[, 1], n = 1)))

# nmi <- sapply(spe_list, function(x) as.numeric(tail(compareClusters(x, func = "NMI")[, 1], n = 1)))

# Visualise pathology annotation and BANKSY cluster on spatial coordinates with the ggspavis package:
library(scater)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(viridis)   
library(ggspavis) # plotSpots


# pal = viridis(n = 16)


# pal = viridis(n = 16)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_bank <- lapply(names(spe_list), function(x) {
  spe_tmp = spe_list[[x]]
  nclu = length(unique(colData(spe_tmp)[,sprintf("%s_smooth", cnm)]))
  plotSpots(spe_tmp, 
            annotate = sprintf("%s_smooth", cnm), 
            size = 0.8, 
            palette = getPalette(nclu),
            in_tissue = NULL) +
    theme(legend.position = "none") +
    labs(title = paste0(x, ",BANKSY clusters"))+ 
    theme(plot.title = element_text(size=33)) + 
    coord_equal()
})

# plot_anno <- lapply(spe_list, function(x) {
#   plotSpots(x, annotate = "clust_annotation", size = 0.8, palette = pal) +
#     theme(legend.position = "none") +
#     labs(title = sprintf("Sample %s: Annotation", x$sample_id[1]))
# })

# plot_list <- c(plot_anno, plot_bank)

pdf("multisample_banksy_lambda1_nucleus.pdf", width = 40, height = 20)
plot_grid(plotlist = plot_bank, ncol = 3, byrow = TRUE)
dev.off()
# plot_bank[[1]]








