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
library(Banksy)
library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)

library(scater)
library(cowplot)
library(ggplot2)



######################################################################### 

xens2 = readRDS(paste0(dir_data_processed, 
                       "xens2.rds")) # spe obj, not too much QC filteration

######################################################################### 
sampleID_unique = as.character(unique(xens2$sample_id))
# > sampleID_unique
# [1] br6197  br5993a br5993b br1735a br1735b br6588 


for(ksample in 1){
  sampleID_tmp = sampleID_unique[ksample]
  se <- xens2[, xens2$sample_id == sampleID_tmp]
  # > se
  # class: SpatialExperiment 
  # dim: 366 94,956  -- use 5% induced points in mNSF
  
  # QC based on total counts
  qcstats <- perCellQCMetrics(se)
  thres <- quantile(qcstats$total, c(0.05, 0.98))
  keep <- (qcstats$total > thres[1]) & (qcstats$total < thres[2])
  se <- se[, keep]
  
  # Normalization to mean library size
  se <- computeLibraryFactors(se)
  aname <- "normcounts"
  assay(se, aname) <- normalizeCounts(se, log = FALSE)
  
  
  # Compute the neighborhood matrices for BANKSY. Setting compute_agf=TRUE computes both the weighted neighborhood mean (
  # ) and the azimuthal Gabor filter (
  # ). The number of spatial neighbors used to compute 
  # and 
  # are k_geom[1]=15 and k_geom[2]=30 respectively. We run BANKSY at lambda=0 corresponding to non-spatial clustering, and lambda=0.2 corresponding to BANKSY for cell-typing.
  
  lambda <- c(0.5, 0.7)
  #lambda: A numeric vector in \in [0,1] specifying a spatial weighting parameter. Larger values incorporate more spatial neighborhood information.
  
  
  k_geom <- c(15, 30)
  #k_geom: An integer scalar specifying the number of neighbors to use. Values \in [15,30] work well.
  
  
  
  se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
  
  
  # Next, run PCA on the BANKSY matrix and perform clustering. Setting use_agf=TRUE uses both 
  # and 
  # to construct the BANKSY matrix.
  set.seed(1000)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
  
  library(Matrix) #downgraded version, to be compatible with Banksy::runBanksyUMAP
  # remotes::install_version("Matrix", version = "1.6-1")
  
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
  se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 1.2)
  
  
  
  # Different clustering runs can be relabeled to minimise their differences with connectClusters:
  se <- Banksy::connectClusters(se)
  #> clust_M1_lam0.2_k50_res1.2 --> clust_M1_lam0_k50_res1.2
  #> 
  # Visualise the clustering output for non-spatial clustering (lambda=0) and BANKSY clustering (lambda=0.2).
  cnames <- colnames(colData(se))
  cnames <- cnames[grep("^clust", cnames)]
  colData(se) <- cbind(colData(se), spatialCoords(se))
  
  saveRDS(se, file = paste0(dir_data_processed_banksy, "se_sample",sampleID_tmp,"_preprosssed_lambda_05_07.rds"))
  
}


saveRDS(sampleID_unique, file=paste0(dir_data_processed_banksy, "sampleID_unique.rds"))

######################################################################### 
### plot the results
######################################################################### 
# plist = list()
sampleID_unique  =readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
pdf("banksy_lambda_05_07.pdf", width = 40, height = 20)
for(ksample in 1:6){
  print(ksample)
  sampleID_tmp = sampleID_unique[ksample]
  se = readRDS(paste0(dir_data_processed_banksy, "se_sample",sampleID_tmp,"_preprosssed_lambda_05_07.rds"))
  cnames <- colnames(colData(se))
  cnames <- cnames[grep("^clust", cnames)]
  colData(se) <- cbind(colData(se), spatialCoords(se))
  
  plot_nsp <- plotColData(se,
                          x = "x_location", y = "y_location",
                          point_size = 0.6, colour_by = cnames[1]
  ) + coord_equal() + ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
  
  plot_bank <- plotColData(se,
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames[2]
  ) + coord_equal()+ ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
  
  # plist[[ksample*2-1]] = plot_nsp
  # plist[[ksample*2]] = plot_bank
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
    
  rdnames <- reducedDimNames(se)
  
  umap_lambdaV1 <- plotReducedDim(se,
                             dimred = grep("UMAP.*lam0.5$", rdnames, value = TRUE),
                             colour_by = cnames[1]
  ) + ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
  
  
  umap_lambdaV2 <- plotReducedDim(se,
                              dimred = grep("UMAP.*lam0.7$", rdnames, value = TRUE),
                              colour_by = cnames[2]
  ) + ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
 
  
  print(plot_grid(plot_nsp , plot_bank , ncol = 2))
  print(plot_grid(
    plot_nsp + facet_wrap(~colour_by),
    plot_bank + facet_wrap(~colour_by),
    ncol = 2
  ))
  
  print(plot_grid(
    umap_lambdaV1,
    umap_lambdaV2,
    ncol = 2
  ))
  
}
dev.off()








