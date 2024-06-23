library(data.table)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(scuttle)
library(Banksy)
library(BiocParallel)
## faster UMAP in uwot 0.2 (used by banksy UMAP)
library(RcppHNSW)

## set wd using rel paths from code dir, where job is launched
setwd("../../")
stopifnot(length(grep(getwd(),pattern="xenium_HYP$",value=T))==1)

## load QCed, genes only, normalized data
hypx <- readRDS("processed-data/04_general_QC_and_normalization/02b-sfe_filtered_genetargonly_log-and-nonlog-norms.RDS")

## set arguments, identical to Yi's arguments used for Banksy cell typing
lam <- 0
res <- 2
gabor <- FALSE
# note that k geom of 6 is suggested in Banksy docs for VISIUM.
# but this is what was used previously so here we'll start
kgeom <- 6
# yi's code used a set seed of 1000 so we'll do that as well to try and hew
# as close as poss to the orig results
monocot <- 1000


## re normalization: Banksy uses  NON-log normcounts (per its documentation, though no
## explanation as to whether log vs nonlog makes an appreciable difference)
## Yi's code used only gene-targeting probes and nonlog norm
## SO: we'll use non-log norm, gene-targeted probes only
## the loaded SFE already contains the normaliezd counts and is subsetted to gene-targeting probes, so nothing to do on these fronts.

## to run Banksy in multisample mode, adjust coordinates so that
## all samples have unique non-overlapping coords (i.e. to prevent
## cells from a given area on two different slides from colliding)
## code from: https://prabhakarlab.github.io/Banksy/articles/batch-correction.html
locs <- spatialCoords(hypx)
locs <- cbind(locs, sample_id = factor(hypx$sample_id))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(hypx)
spatialCoords(hypx) <- locs
##

## save SFE with normed counts and adjusted coords
saveRDS(hypx,"processed-data/05_Banksy_M0lam0_res2_multisamp/01-sfe-in_staggered-coords_genetarg-only_NONlog-norm.RDS")

# now specify which Assays() slot to use for banksy calls
useassay <- "normcounts"

# aaand run Banksy
hypx <- computeBanksy(hypx,
                      assay_name=useassay,
                      compute_agf = gabor,
                      k_geom = kgeom,
                      seed=monocot)

## run Banksy PCA for subseequent input to clustering
## CRITICAL HERE: specify "group" as the sample_id -- i.e., each
## sample's Banksy PCA will be computed separately, making this
## equivalent to running Banksy on a list of SFEs with each entry as
## one sample, the other approach their vignettes uses
hypx <- runBanksyPCA(hypx,
                         assay_name=useassay,
                         use_agf=gabor,
                         lambda=lam,
                         group="sample_id",
                         seed=monocot)

## run Banksy clustering -- Yi used the default, Leiden clustering.
## this will take a longgggg time.
hypx <- clusterBanksy(hypx,
                      assay_name=useassay,
                      use_agf = gabor,
                      lambda = lam,
                      resolution = res,
                      group="sample_id",
                      seed = monocot)


## make the clustering result non-numeric for ease of plotting later;
## extract the clustering result, rather than saving a whole nother SFE
colData(hypx)$clust_M0_lam0_k50_res2 <- as.character(colData(hypx)$clust_M0_lam0_k50_res2)
colData(hypx)$clust_M0_lam0_k50_res2 <- paste0("X",colData(hypx)$clust_M0_lam0_k50_res2)

m0lam0res2 <- as.data.table(colData(hypx),keep.rownames=T)

m0lam0res2 <- m0lam0res2[,.(rn,clust_M0_lam0_k50_res2)]
setnames(m0lam0res2,c("rn","M0lam0_leiden_multisamp_res2"))
fwrite(m0lam0res2,
       "processed-data/05_Banksy_M0lam0_res2_multisamp/01-BanksyClusts_M0lam0_leiden_multisamp.txt",
       sep='\t',
       quote=F)

## reproducibility info
sessionInfo()
