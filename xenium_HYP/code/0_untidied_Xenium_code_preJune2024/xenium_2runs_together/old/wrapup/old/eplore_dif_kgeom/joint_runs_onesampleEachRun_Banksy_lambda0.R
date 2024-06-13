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

# dim(spe_lists_old[[1]])
# dim(spe_list_new[[1]])
# dim(spe_list_new[[2]])
# dim(spe_list_new[[3]])

# spe_list = c(spe_lists_old, spe_list_new)
genes_shared = intersect(rownames(spe_lists_old[[1]]), rownames(spe_list_new[[1]]))
# > length(genes_shared)
# [1] 366

for(k in 1:length(spe_lists_old)){
  spe_lists_old[[k]] = spe_lists_old[[k]][genes_shared,]
}

for(k in 1:length(spe_list_new)){
  spe_list_new[[k]] = spe_list_new[[k]][genes_shared,]
}
# spe_list = c(list(spe_lists_old[[1]][genes_shared,]),
# list(spe_list_new[[1]][genes_shared,]))

colData(spe_lists_old[[1]]) = colData(spe_lists_old[[1]])[,colnames(colData(spe_list_new[[1]]))]
spe_list = list(spe_lists_old[[1]], spe_list_new[[1]])

###########
compute_agf <- FALSE
# k_geom=6
k_geom=12

# for(k_geom in c(6,12)){
  print("k_geom")
  print(k_geom)
  
  spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)
  
  spe_joint <- do.call(cbind, spe_list)

  
  # When running multi-sample BANKSY PCA, the group argument may be provided. This specifies the grouping variable for the cells or spots across the samples. Features belonging to cells or spots corresponding to each level of the grouping variable will be z-scaled separately. In this case, sample_id in colData(spe_joint) gives the grouping based on the sample of origin.
  lambda <- 0
  use_agf <- FALSE
  spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)
  
  # Run UMAP on the BANKSY embedding:
  spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
  saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda0_oneSampleEachRun.rds"))
  
  # Finally, obtain cluster labels for spots across all 4 samples.
  res <- 0.6
  spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
  cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
  saveRDS(spe_joint, file = paste0("spe_joint_banksy_lambda0_oneSampleEachRun.rds"))
  
  
  # Once joint clustering is performed, we split the samples into their own SpatialExperiment objects:
  sample_names= c(sample_names_old, sample_names_new)
  spe_list <- lapply(sample_names, function(x) spe_joint[, spe_joint$sample_id == x])

  # 
  # As an optional step, we smooth the cluster labels of each sample separately. This can be done if smooth spatial domains are expected in the biological sample or tissue in question.
  spe_list <- lapply(spe_list, smoothLabels, cluster_names = cnm, k = 6L, verbose = FALSE)
  names(spe_list) <- paste0("sample_", sample_names)
  
  saveRDS(spe_list, file = paste0("spe_list_banksy_lambda0_oneSampleEachRun0.rds"))
  

  
  ######################################################################### 
  ###################################################################### 
  ## do the plots
  ######################################################################### 
  se = spe_joint
  # unique(colData(se)[,"sample_id"])
  colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
  unique(colData(se)[,"sample_id"] )
  
  
  plt <- function(se_, 
                  grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                  grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
    
    cnames <- colnames(colData(se_))
    
    cnames <- cnames[grep("^clust", cnames)]
    print(cnames)
    
    colData(se_) <- cbind(colData(se_), spatialCoords(se_))
    
    plot_bank1 <- plotColData(se_[,grp1],
                              x = "x_location", y = "y_location",
                              point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
    ) + coord_equal() + 
      ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
      theme(plot.title = element_text(size=33)) + 
      facet_wrap(~ shape_by, ncol = 3)+   
      scale_shape_manual(values=seq(0,15))
    
    plot_bank2 <- plotColData(se_[,grp2],
                              x = "x_location", y = "y_location",
                              point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
    ) + coord_equal() + 
      ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
      theme(plot.title = element_text(size=33)) + 
      facet_wrap(~ shape_by, ncol = 3)+   
      scale_shape_manual(values=seq(0,15))
    
    
    plot_bank_byClu1 <- plotColData(se_[,grp1],
                                    x = "x_location", y = "y_location",
                                    point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
    ) + coord_equal()+ 
      ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
      theme(plot.title = element_text(size=33)) + 
      facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
      scale_shape_manual(values=seq(0,15))
    
    plot_bank_byClu2 <- plotColData(se_[,grp2],
                                    x = "x_location", y = "y_location",
                                    point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
    ) + coord_equal()+ 
      ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
      theme(plot.title = element_text(size=33)) + 
      facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
      scale_shape_manual(values=seq(0,15))
    
    
    # Visualize UMAPs of the non-spatial and BANKSY embedding:
    rdnames <- reducedDimNames(se_)
    
    umap_1 <- plotReducedDim(se_[,grp1],
                             dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                             colour_by = cnames, shape_by = "sample_id"
    ) + 
      ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
      theme(plot.title = element_text(size=33)) + 
      facet_wrap(~shape_by, ncol = 2)+   
      scale_shape_manual(values=seq(0,15))
    
    umap_2<- plotReducedDim(se_[,grp2],
                            dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                            colour_by = cnames, shape_by = "sample_id"
    ) + 
      ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
      theme(plot.title = element_text(size=33)) + 
      facet_wrap(~shape_by, ncol = 2)+   
      scale_shape_manual(values=seq(0,15))
    
    print(plot_bank1)  
    print(plot_bank2)
    
    print(plot_bank_byClu1)
    print(plot_bank_byClu2)
    
    print(umap_1)
    print(umap_2)
    
  }
  library(scater)
  pdf(paste0("multisample_banksy_2runs_oneSampleEachRun_lambda_0.pdf"), width = 40, height = 100/4)
  print(plt(se[,]))
  # plt(se[,])
  dev.off()
  
# }
# 
# 

