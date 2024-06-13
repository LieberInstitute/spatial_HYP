# ml conda_R/4.3
# R

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"

dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/wrapup/one_sample_banksy/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/wrapup/one_sample_banksy/"

dir.create("/dcs04/hansen/data/ywang/xenium/wrapup/")
dir.create("/dcs04/hansen/data/ywang/xenium/wrapup/one_sample_banksy")

setwd(dir_out)

###########
library(Seurat)
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

ksample=1

# for(ksample in 1){
  sampleID_tmp = sampleID_unique[ksample]
  se <- xens2[, xens2$sample_id == sampleID_tmp]
 
  # QC based on total counts
  qcstats <- perCellQCMetrics(se)
  thres <- quantile(qcstats$total, c(0.05, 0.98))
  keep <- (qcstats$total > thres[1]) & (qcstats$total < thres[2])
  se <- se[, keep]
  
  # Normalization to mean library size
  se <- computeLibraryFactors(se)
  aname <- "normcounts"
  assay(se, aname) <- normalizeCounts(se, log = FALSE)
  
  lambda <- c(1)
  
  k_geom <- c(30)
  
  se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
  
  # Next, run PCA on the BANKSY matrix and perform clustering. Setting use_agf=TRUE uses both 
  # and 
  # to construct the BANKSY matrix.
  set.seed(1000)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
  
  library(Matrix) #downgraded version, to be compatible with Banksy::runBanksyUMAP
  # remotes::install_version("Matrix", version = "1.6-1")
  
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
  
  saveRDS(se, file = paste0(dir_data_processed_banksy, 
                                  "se_sample",sampleID_tmp,
                            "_preprosssed_lambda1_UMAP.rds"))
  
  
  