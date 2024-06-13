lambda=1
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
spe_list_fullData = readRDS(paste0("spe_list_banksy_lambda1.rds"))

########### keep ARC-domain spots
spe_list = list()
for(sample in names(spe_list_fullData)){
  
  cellClu = colData(spe_list_fullData[[sample]][,])[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
  # se_seurat[,cellClu=="4"]
  spe_list[[sample]] = spe_list_fullData[[sample]][,cellClu=="3"]
  
}
saveRDS(spe_list, file="spe_list_banksy_lambda1_arcOnly.rds")

###########
spe_list  = readRDS("spe_list_banksy_lambda1_arcOnly.rds")

###########
compute_agf <- FALSE
k_geom <- 6
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname ,compute_agf = compute_agf, k_geom = k_geom)

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
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda1_arcOnly.rds"))
spe_list  = readRDS("spe_list_banksy_lambda1_arcOnly.rds")

#######
spe_joint = readRDS("spe_joint_banksy_lambda1_arcOnly.rds")


# Finally, obtain cluster labels for spots across all 4 samples.
res <- 0.2
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda1_arcOnly_res02.rds"))

unique(cnm)
unique(colData(spe_joint)[cnm])
# clust_M0_lam1_k50_res0.2
# <factor>
#   br6197.aaaaiadp-1                        1
# br6197.aaahimpg-1                        2
# br6197.aabmlfna-1                        4
# br6197.acnekieb-1                        3
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

saveRDS(spe_list, file = paste0("spe_list_banksy_lambda1_arcOnly_res02.rds"))


######################################################################### 
# Parsing BANKSY output
######################################################################### 
spe_list = readRDS("spe_list_banksy_lambda1_arcOnly_res02.rds")
use_agf <- FALSE
lambda <- 1
res <- 0.2
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
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_bank <- lapply(names(spe_list), function(x) {
  spe_tmp = spe_list[[x]]
  
  colData(spe_tmp)[,"clust_M0_lam1_k50_res0.6"]=NULL
  nclu = length(unique(colData(spe_tmp)[,cnm]))
  
  # nclu = length(unique(colData(spe_tmp)[,sprintf("%s_smooth", cnm)]))
  plotSpots(spe_tmp, 
            annotate = (cnm), 
            
            # annotate = sprintf("%s_smooth", cnm), 
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

pdf("multisample_banksy_lambda1_arcOnly_res02.pdf", width = 40, height = 20)
plot_grid(plotlist = plot_bank, ncol = 3, byrow = TRUE)
dev.off()
# # plot_bank[[1]]
# 







