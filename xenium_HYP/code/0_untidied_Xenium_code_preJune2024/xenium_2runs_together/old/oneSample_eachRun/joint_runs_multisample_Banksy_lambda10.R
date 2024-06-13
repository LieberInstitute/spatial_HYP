######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# setwd(dir_out)
###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
aname <- "normcounts"

###########
sample_names_old = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))

# To run BANKSY across multiple samples, we first compute the BANKSY neighborhood feature matrices for each sample separately. We use k_geom=6 corresponding to the first-order neighbors in 10x Visium assays (k_geom=18 corresponding to first and second-order neighbors may also be used).
spe_lists_old = readRDS(paste0(dir_data_processed_banksy,"spe_list.rds"))
spe_list_new = readRDS(paste0(dir_data_processed_banksy_new,"spe_list.rds"))

dim(spe_lists_old[[1]])
dim(spe_list_new[[1]])

genes_shared = intersect(rownames(spe_lists_old[[1]]), rownames(spe_list_new[[1]]))
# > length(genes_shared)
# [1] 366

spe_list = c(list(spe_lists_old[[1]][genes_shared,]),
             list(spe_list_new[[1]][genes_shared,]))

###########
compute_agf <- FALSE
k_geom <- 6
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)

# We then merge the samples to perform joint dimensional reduction and clustering:
spe_joint <- do.call(cbind, spe_list)
# rm(spe_list)
# gc()


# When running multi-sample BANKSY PCA, the group argument may be provided. This specifies the grouping variable for the cells or spots across the samples. Features belonging to cells or spots corresponding to each level of the grouping variable will be z-scaled separately. In this case, sample_id in colData(spe_joint) gives the grouping based on the sample of origin.
lambda <- 0
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)

# Run UMAP on the BANKSY embedding:
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda0.rds"))

# Finally, obtain cluster labels for spots across all 4 samples.
res <- 0.6
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda0.rds"))

# To compare clusters visually in the next section, run connectClusters:
# spe_joint <- connectClusters(spe_joint) # skipped this line, due to error messages 
# connectClusters: Relabel cluster labels across parameter runs to maximise their similarity.
# NA --> clust_M0_lam0.2_k50_res0.6
# Error in `[.data.frame`(clust_df, , curr_many) : 
  # undefined columns selected


# Once joint clustering is performed, we split the samples into their own SpatialExperiment objects:
sample_names= c(sample_names_old[1], sample_names_new[1])
spe_list <- lapply(sample_names, function(x) spe_joint[, spe_joint$sample_id == x])
# rm(spe_joint)
# gc()
# 
# As an optional step, we smooth the cluster labels of each sample separately. This can be done if smooth spatial domains are expected in the biological sample or tissue in question.
spe_list <- lapply(spe_list, smoothLabels, cluster_names = cnm, k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_names)

saveRDS(spe_list, file = paste0("spe_list_banksy_lambda0.rds"))

######################################################################### 
# Parsing BANKSY output
######################################################################### 
spe_list = readRDS("spe_list_banksy_lambda0.rds")
use_agf <- FALSE
lambda <- 0
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

pdf("multisample_banksy_lambda0.pdf", width = 40, height = 20)
plot_grid(plotlist = plot_bank, ncol = 2, byrow = TRUE)
dev.off()
# plot_bank[[1]]








