######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R
###########
res <- 4
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
# dir_data_processed =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
dir_data_processed_banksy =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# setwd(dir_out)
###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"
# dir.create(dir_data_processed_banksy)
# dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)
###########

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
aname <- "normcounts"

###########
dir_allSamples_old = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20231019__161052__101923_KMon101923/"
list_sample_old = list.files(dir_allSamples_old)

dir_allSamples_new= "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20240223__174743__022324_KMon120823/"
list_sample_new = list.files(dir_allSamples_new)

###########
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")

# To run BANKSY across multiple samples, we first compute the BANKSY neighborhood feature matrices for each sample separately. We use k_geom=6 corresponding to the first-order neighbors in 10x Visium assays (k_geom=18 corresponding to first and second-order neighbors may also be used).
spe_lists_old = readRDS(paste0(dir_data_processed_banksy,"spe_list_kp.rds"))
spe_list_new = readRDS(paste0(dir_data_processed_banksy_new,"spe_list_kp.rds"))

dim(spe_lists_old[[1]])
dim(spe_list_new[[1]])

genes_shared = intersect(rownames(spe_lists_old[[1]]), rownames(spe_list_new[[1]]))
# > length(genes_shared)
# [1] 366
genes_shared <- grep(genes_shared,pattern="NegControlProbe|NegControlCodeword|DeprecatedCodeword|BLANK",value = T,invert = T)

spe_list = c(spe_lists_old, spe_list_new)

for(k in 1:length(spe_list)){
  spe_list[[k]] =  spe_list[[k]][genes_shared,]
}


print("00")

###########
compute_agf <- FALSE
k_geom <- 6
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)

# We then merge the samples to perform joint dimensional reduction and clustering:
spe_joint <- do.call(cbind, spe_list)
# rm(spe_list)
# gc()

print("11")

# When running multi-sample BANKSY PCA, the group argument may be provided. This specifies the grouping variable for the cells or spots across the samples. Features belonging to cells or spots corresponding to each level of the grouping variable will be z-scaled separately. In this case, sample_id in colData(spe_joint) gives the grouping based on the sample of origin.
lambda <- 0
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)

# Run UMAP on the BANKSY embedding:
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
saveRDS(spe_joint, file = paste0("spe_joint_2runs_banksy_lambda0_res", res, ".rds"))
print("22")

# Finally, obtain cluster labels for spots across all 4 samples.


spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
saveRDS(spe_joint, file = paste0("spe_joint_2runs_banksy_lambda0_res", res, ".rds"))
print("33")

# To compare clusters visually in the next section, run connectClusters:
# spe_joint <- connectClusters(spe_joint) # skipped this line, due to error messages 
# connectClusters: Relabel cluster labels across parameter runs to maximise their similarity.
# NA --> clust_M0_lam0.2_k50_res0.6
# Error in `[.data.frame`(clust_df, , curr_many) : 
  # undefined columns selected



#############
#############
spe_joint = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", res, ".rds"))


# Once joint clustering is performed, we split the samples into their own SpatialExperiment objects:
sample_names= c(sample_names_old, sample_names_new)
spe_list <- lapply(sample_names, function(x) spe_joint[, spe_joint$sample_id == x])
# rm(spe_joint)
# gc()

# As an optional step, we smooth the cluster labels of each sample separately. This can be done if smooth spatial domains are expected in the biological sample or tissue in question.
spe_list <- lapply(spe_list, smoothLabels, cluster_names = cnm, k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_names)

saveRDS(spe_list, file = paste0("spe_list_2runs_banksy_lambda0_res", res, ".rds"))
print("44")

######################################################################### 
# Parsing BANKSY output
######################################################################### 
spe_list = readRDS("spe_list_2runs_banksy_lambda0_res", res, ".rds")
use_agf <- FALSE
lambda <- 0
# res <- 0.6
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
print("55")

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
coul <- brewer.pal(9, "Set1")
# Add more colors to this palette :
# coul <- colorRampPalette(coul)(nclu)
nclu = length(unique(colData(spe_joint)[, cnm]))

plot_bank <- lapply(names(spe_list), function(x) {
  spe_tmp = spe_list[[x]]
  # nclu = length(unique(colData(spe_tmp)[,sprintf("%s_smooth", cnm)]))
  # nclu=length(unique(colData(spe_tmp)[,cnm]))
  # display.brewer.all(n=10, exact.n=FALSE)

  
  # getPalette = colorRampPalette(brewer.pal(nclu, "greenArmytage"))
  getPalette = colorRampPalette(coul)
  
  plotSpots(spe_tmp, 
            annotate = cnm, 
            
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

pdf("multisample_2runs_banksy_lambda0_res", res, ".pdf", width = 60, height = 100)
plot_grid(plotlist = plot_bank, ncol = 3, byrow = TRUE)
dev.off()
# plot_bank[[1]]


# print("66")






