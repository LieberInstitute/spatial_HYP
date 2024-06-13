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
                       "xens2.rds")) # spe obj, not too much QC filteration
sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))


harmo = readRDS( paste0(dir_data_processed, 
                        "xens2_harmony.rds"))


## seurat louvain multilevel refinement clustering? 
# xens3  <- as.Seurat(harmo) # added by Yi
# assay(spe, "counts") <- NULL

# assays(spe)$counts  = Embeddings(xens3, reduction = "HARMONY")
# Embeddings(xens3)

# dim(spe)
# dim(xens3)
# [1]    366 407710
# [1]    366 407710
# colnames(xens3)


spe_harmonyEmb <- SpatialExperiment(
  assay = list(counts = t(Embeddings(xens3, reduction = "HARMONY"))), 
  colData = colData(harmo), 
  spatialCoordsNames = c("x_position", "y_position"))

saveRDS(spe_harmonyEmb, file = "spe_harmonyEmb.rds")

######################################################################### 
#################### step1 create spe_list
######################################################################### 
imgData(spe_harmonyEmb) <- NULL
assay(spe_harmonyEmb, "logcounts") <- NULL
reducedDims(spe_harmonyEmb) <- NULL
rowData(spe_harmonyEmb) <- NULL
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
spe_list <- lapply(sample_names, function(x) spe_harmonyEmb[, spe_harmonyEmb$sample_id == x])
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


saveRDS(spe_list, file = paste0(dir_data_processed_banksy,"spe_list_harmony.rds"))



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
sample_names  =readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))

# To run BANKSY across multiple samples, we first compute the BANKSY neighborhood feature matrices for each sample separately. We use k_geom=6 corresponding to the first-order neighbors in 10x Visium assays (k_geom=18 corresponding to first and second-order neighbors may also be used).
spe_list = readRDS(paste0(dir_data_processed_banksy,"spe_list_harmony.rds"))

###########
compute_agf <- FALSE
k_geom <- 6
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)

# We then merge the samples to perform joint dimensional reduction and clustering:
spe_joint <- do.call(cbind, spe_list)
rm(spe_list)
gc()


# When running multi-sample BANKSY PCA, the group argument may be provided. This specifies the grouping variable for the cells or spots across the samples. Features belonging to cells or spots corresponding to each level of the grouping variable will be z-scaled separately. In this case, sample_id in colData(spe_joint) gives the grouping based on the sample of origin.
lambda <- 0.2
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda02_onHarmonyEmbeddings.rds"))


# spe_joint = readRDS("spe_joint_banksy_lambda02_onHarmonyEmbeddings.rds")

# Run UMAP on the BANKSY embedding:
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda02_onHarmonyEmbeddings.rds"))

# spe_joint = readRDS("spe_joint_banksy_lambda02_onHarmonyEmbeddings.rds")

# Finally, obtain cluster labels for spots across all 4 samples.
res <- 0.6
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda02_onHarmonyEmbeddings.rds"))

# To compare clusters visually in the next section, run connectClusters:
# spe_joint <- connectClusters(spe_joint)

# Once joint clustering is performed, we split the samples into their own SpatialExperiment objects:
spe_list <- lapply(sample_names, function(x) spe_joint[, spe_joint$sample_id == x])
# rm(spe_joint)
# gc()

# As an optional step, we smooth the cluster labels of each sample separately. 
# This can be done if smooth spatial domains are expected in the biological sample or tissue in question.
spe_list <- lapply(spe_list, smoothLabels, cluster_names = cnm, k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_names)

saveRDS(spe_list, file = paste0("spe_list_banksy_lambda02_onHarmonyEmbeddings.rds"))


######################################################################### 
# Parsing BANKSY output
######################################################################### 
# We can compare BANKSY clusters to pathology annotations using several cluster comparison measures such as the adjusted Rand index (ARI) or normalized mutual information (NMI) with compareClusters.
# 
# ari <- sapply(spe_list, function(x) as.numeric(tail(compareClusters(x, func = "ARI")[, 1], n = 1)))

# nmi <- sapply(spe_list, function(x) as.numeric(tail(compareClusters(x, func = "NMI")[, 1], n = 1)))

# Visualise pathology annotation and BANKSY cluster on spatial coordinates with the ggspavis package:
  
# Use scater:::.get_palette('tableau10medium')
# pal <- c(
#     "#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9",
#              "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA"
# )

# colourCount = length(unique(mtcars$hp))

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

# plot_bank <- lapply(names(spe_list), function(x) {
#   plotSpots(spe_list[[x]], annotate = sprintf("%s_smooth", cnm), size = 0.8, palette = getPalette(16),
#             in_tissue = NULL) +
#     theme(legend.position = "none") +
#     labs(title = paste0(x, "BANKSY clusters"))
# })

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

# plot_grid(plotlist = plot_bank, ncol = 3, byrow = TRUE)

pdf("multisample_banksy_onHarmonyEmbeddings_lambda02.pdf", width = 40, height = 20)
plot_grid(plotlist = plot_bank, ncol = 3, byrow = TRUE)
dev.off()
























######################################################################### 
# sampleID_unique = as.character(unique(xens2$sample_id))
# > sampleID_unique
# [1] br6197  br5993a br5993b br1735a br1735b br6588 

for(ksample in 6){
  sampleID_tmp = sampleID_unique[ksample]

  se = readRDS(paste0(dir_data_processed_banksy, "se_sample",sampleID_tmp,"_preprosssed_lambda_099_1.rds"))
  # colData(se)
  # se_seurat = readRDS( paste0(dir_data_processed_banksy, "se_seurat_sample",sampleID_tmp,"_markerGenesBanksy.rds"))
  # saveRDS(se_seurat, file = paste0(dir_data_processed_banksy, "se_seurat_sample",sampleID_tmp,"_markerGenesBanksy.rds"))
  harmo <- harmony::RunHarmony(se, 
                               lambda=NULL, 
                               group.by.vars=c("brnum","slide"), 
                               ncores=1)
  
}





######################################################################### 
### plot the results
######################################################################### 
