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
k_geom <- 30
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname ,
                   compute_agf = compute_agf, k_geom = k_geom)

# We then merge the samples to perform joint dimensional reduction and clustering:
spe_joint <- do.call(cbind, spe_list)
# rm(spe_list)
# gc()


# When running multi-sample BANKSY PCA, the group argument may be provided. This specifies the grouping variable for the cells or spots across the samples. Features belonging to cells or spots corresponding to each level of the grouping variable will be z-scaled separately. In this case, sample_id in colData(spe_joint) gives the grouping based on the sample of origin.
lambda <- 1
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, 
                          group = "sample_id", seed = 1000)

# Run UMAP on the BANKSY embedding:
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda1_arcOnly_30neighbours.rds"))
# spe_list  = readRDS("spe_list_banksy_lambda1_arcOnly.rds")

#######
spe_joint = readRDS("spe_joint_banksy_lambda1_arcOnly_30neighbours.rds")


# Finally, obtain cluster labels for spots across all 4 samples.
res <- 0.4
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, 
                           resolution = res, seed = 1000)
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda1_arcOnly_res04_30neighbours.rds"))

spe_joint = readRDS("spe_joint_banksy_lambda1_arcOnly_res04_30neighbours.rds")
# colData(spe_joint)[cnm]=NULL
# colData(spe_joint)[cnm]=NULL
# colData(spe_joint)[cnm]=NULL

unique(cnm)
# 14 clusters
unique(colData(spe_joint)[cnm])


# Once joint clustering is performed, we split the samples into their own SpatialExperiment objects:
spe_list <- lapply(sample_names, function(x) spe_joint[, spe_joint$sample_id == x])
# rm(spe_joint)
# gc()
# 
# As an optional step, we smooth the cluster labels of each sample separately. This can be done if smooth spatial domains are expected in the biological sample or tissue in question.
spe_list <- lapply(spe_list, smoothLabels, cluster_names = cnm, k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_names)

saveRDS(spe_list, file = paste0("spe_list_banksy_lambda1_arcOnly_res04_30neighbours.rds"))

######################################################################### 
# Parsing BANKSY output
######################################################################### 
spe_list = readRDS("spe_list_banksy_lambda1_arcOnly_res04_30neighbours.rds")
use_agf <- FALSE
lambda <- 1
res <- 0.4
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


########### load multi-sample Banksy result
# spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))

spe_joint  = readRDS(paste0("spe_joint_banksy_lambda1_arcOnly_res04_30neighbours.rds"))
# spe_joint_banksy_lambda1_res2.rds
se = spe_joint
# Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

head(colData(se))

colData(se)$clust_M0_lam1_k50_res0.6 = NULL
colData(se)$clust_M0_lam1_k50_res0.6_smooth = NULL

######################################################################### 
## do the plots
######################################################################### 
plt <- function(){
  
  cnames <- colnames(colData(se))
  

  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
  colData(se) <- cbind(colData(se), spatialCoords(se))
  
  plot_bank <- plotColData(se,
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)
  
  plot_bank_byClu <- plotColData(se,
                                 x = "x_location", y = "y_location",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = 6)
  
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
  rdnames <- reducedDimNames(se)
  
  umap_ <- plotReducedDim(se,
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 3)
  
  print(plot_bank)
  print(plot_bank_byClu)
  
  print(umap_)
}
lambdaName=1
# sampleID_unique  = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
pdf(paste0("multisample_banksy_lambda_", lambdaName, "_arcOnly_res04_30neighbours.pdf"), width = 40, height = 60)
plt()
dev.off()



pdf(paste0("multisample_banksy_lambda_", lambdaName, "_largerPage_arcOnly_res04_30neighbours.pdf"), width = 80, height = 60)
plt()
dev.off()

pdf(paste0("multisample_banksy_lambda_", lambdaName, "_smallerPage_arcOnly_res04_30neighbours.pdf"), width = 20, height = 30)
plt()
dev.off()
# 
# 
# ######################################################################### 
# # Parsing BANKSY output
# ######################################################################### 
# spe_list = readRDS("spe_list_banksy_lambda1_arcOnly_30neighbours_res04.rds")
# use_agf <- FALSE
# lambda <- 1
# res <- 0.4
# cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
# 
# # We can compare BANKSY clusters to pathology annotations using several cluster comparison measures such as the adjusted Rand index (ARI) or normalized mutual information (NMI) with compareClusters.
# # 
# # ari <- sapply(spe_list, function(x) as.numeric(tail(compareClusters(x, func = "ARI")[, 1], n = 1)))
# 
# # nmi <- sapply(spe_list, function(x) as.numeric(tail(compareClusters(x, func = "NMI")[, 1], n = 1)))
# 
# # Visualise pathology annotation and BANKSY cluster on spatial coordinates with the ggspavis package:
# library(scater)
# library(ggpubr)
# library(cowplot)
# library(RColorBrewer)
# library(viridis)   
# library(ggspavis) # plotSpots
# 
# 
# # pal = viridis(n = 16)
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# 
# plot_bank <- lapply(names(spe_list), function(x) {
#   spe_tmp = spe_list[[x]]
#   colData(spe_tmp)[,"clust_M0_lam1_k50_res0.6"]=NULL
#   # colData(spe_tmp)[,"clust_M0_lam1_k50_res0.6"]=NULL
#   nclu = length(unique(colData(spe_tmp)[,cnm]))
#   
#   # nclu = length(unique(colData(spe_tmp)[,sprintf("%s_smooth", cnm)]))
#   plotSpots(spe_tmp, 
#             annotate = (cnm), 
#             
#             # annotate = sprintf("%s_smooth", cnm), 
#             size = 0.8, 
#             palette = getPalette(nclu),
#             in_tissue = NULL) +
#     theme(legend.position = "none") +
#     labs(title = paste0(x, ",BANKSY clusters"))+ 
#     theme(plot.title = element_text(size=33)) + 
#     coord_equal()
# })
# 
# # plot_anno <- lapply(spe_list, function(x) {
# #   plotSpots(x, annotate = "clust_annotation", size = 0.8, palette = pal) +
# #     theme(legend.position = "none") +
# #     labs(title = sprintf("Sample %s: Annotation", x$sample_id[1]))
# # })
# 
# # plot_list <- c(plot_anno, plot_bank)
# 
# pdf("multisample_banksy_lambda1_arcOnly_res04_30neighbours.pdf", width = 40, height = 20)
# plot_grid(plotlist = plot_bank, ncol = 3, byrow = TRUE)
# dev.off()
# # # plot_bank[[1]]
# # 
# 
# 
# 
# 
# 
# 
# 
