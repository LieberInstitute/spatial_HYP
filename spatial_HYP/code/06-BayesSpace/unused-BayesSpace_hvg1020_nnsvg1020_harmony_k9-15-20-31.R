library(data.table)
library(Biostrings)
library(SpatialExperiment)
library(scater) # addPerCellQC
library(scran)
library(BayesSpace)

#### Load SPE containing Harmony dimensionality reductions for 10%ile and 20%ile HVGs, 10%ile and 20%ile mean-rank nnSVGs, each with bayesspace k=2.
setwd("/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/")
hyp <- readRDS("data/04-feature_selection/02-hypfiltered_4featureset-pca-umap-harmonydefault-harmonylambdanull-mnn30.RDS")
## make certain colnames are the same as $key or we won't have unique spot IDs to pair the assignments back to!
if(sum(colnames(hyp)==hyp$key)!=ncol(hyp)){colnames(hyp) <- hyp$key}

### build table with job-wise variables (names of dimensionality reductions to use in bayesspace and value to pass to q in BayesSpace from among values of 9/15/20/31)

p <- c(rep(c("HARMONYdflt_hvg10","HARMONYdflt_hvg20","HARMONYdflt_nnsvg10","HARMONYdflt_nnsvg20","HARMONYlmbna_hvg10","HARMONYlmbna_hvg20","HARMONYlmbna_nnsvg10","HARMONYlmbna_nnsvg20","mnn30_hvg10","mnn30_hvg20","mnn30_nnsvg10","mnn30_nnsvg20"),4))
q <- c(rep(9,12),rep(15,12),rep(20,12),rep(31,12))
pq <- cbind(p,q)
# so that rownames still end up == to row index after unique'ing:
pq <- as.data.table(unique(pq))
### ^ 48 runs (4 q values per feature-reduction x 4 feature sets x 3 reductions (harmony defaults, harmony lambda=null, mnn30))
pq[,q:=as.numeric(q)]
pq <- as.data.frame(pq)
stopifnot(nrow(unique(pq))==nrow(pq))
rm(p,q)

### pull params for task number
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
curdimred <- pq[i,1]
curq <- pq[i,2]

## important from slack: need to run the next two lines or will get "subscript contains invalid names"

colData(hyp)$row <- colData(hyp)$array_row
colData(hyp)$col <- colData(hyp)$array_col


### according to https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/03_BayesSpace/01_BayesSpace.R , don't use spatial preprocess. in order to do this you have to reset metadata.
### HOWEVER, BayesSpace now throws an error if you don't actually run spatialpreproc. also, the documentation says that it ADDS metadata, which shouldn't be equiv to overwriting it. The "metadata" object in the HYP SPE is an empty list, so nothing will be modified. (Maybe LIBD used to place something in the metadata entry for something but evidently not now).

spatialPreprocess(hyp,platform="Visium",skip.PCA=T)


# Main call
# defaults to 50k reps, LIBD uses 10k

tmp <- spatialCluster(sce = hyp,init.method = "mclust",use.dimred = curdimred,q = curq,platform = "Visium",nrep=10000)

# pull out clusters to save in much tinier file than an SPE
tmp2 <- as.data.table(colData(tmp)[,c("key","spatial.cluster")])
setnames(tmp2,c("rn",paste0("BShar_",curdimred,"_k",curq)))

fwrite(tmp2,paste0("data/06-BayesSpace/01-bayesspace_k9-15-20-31_out/BSpace_k",curq,"_",curdimred,".txt"),sep='\t',quote=F)
rm(list=ls())
gc(full=T)

library(sessioninfo)
sessionInfo()
session_info()
