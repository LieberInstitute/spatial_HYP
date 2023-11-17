library(data.table)
library(Biostrings)
library(SpatialExperiment)
library(scater) # addPerCellQC
library(scran)
library(BayesSpace)

#### Load SPE containing Harmony dimensionality reductions for 10%ile and 20%ile HVGs, 10%ile and 20%ile mean-rank nnSVGs, each with bayesspace k=2.

hyp <- readRDS("../../data/procesing/hypfiltered_pca-umap-harmonydefault-harmonylambdanull.RDS")
## make certain colnames are the same as $key or we won't have unique spot IDs to pair the assignments back to!
colnames(hyp) <- hyp$key

### names of dimensionality reductions to use in bayesspace

p <- c("HARMONYdflt_hvg10","HARMONYdflt_hvg20","HARMONYdflt_nnsvg10","HARMONYdflt_nnsvg20","HARMONYlmbna_hvg10","HARMONYlmbna_hvg20","HARMONYlmbna_nnsvg10","HARMONYlmbna_nnsvg20","mnn30_hvg10","mnn30_hvg20","mnn_nnsvg10","mnn_nnsvg20")

### pull params for task number

i <- as.numeric(Sys.getenv("SGE_TASK_ID"))
curdimred <- p[i]

## important from slack: need to run the next two lines or will get "subscript contains invalid names"

colData(hyp)$row <- colData(hyp)$array_row
colData(hyp)$col <- colData(hyp)$array_col


### according to https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/03_BayesSpace/01_BayesSpace.R , don't use spatial preprocess. in order to do this you have to reset metadata.
### HOWEVER, BayesSpace now throws an error if you don't actually run spatialpreproc. also, the documentation says that it ADDS metadata, which shouldn't be equiv to overwriting it. The "metadata" object in the HYP SPE is an empty list, so nothing will be modified. (Maybe LIBD used to place something in the metadata entry for something but evidently not now).

spatialPreprocess(hyp,platform="Visium",skip.PCA=T)


# Main call
# defaults to 50k reps, LIBD uses 10k, we'll use 20k here

tmp <- spatialCluster(sce = hyp,init.method = "mclust",use.dimred = curdimred,q = 2,platform = "Visium",nrep=20000)

# pull out clusters to save in much tinier file than an SPE
tmp2 <- as.data.table(cbind(rownames(colData(tmp)),colData(tmp)$spatial.cluster))
setnames(tmp2,c("rn",paste0("BShar_",curname,"_",curdimred,"_k2")))

fwrite(tmp2,paste0("harmony_bs_k2_out/BSpace_k2_",curdimred,".txt"),sep='\t',quote=F)
rm(list=ls())
gc(full=T)

library(sessioninfo)
session_info()
