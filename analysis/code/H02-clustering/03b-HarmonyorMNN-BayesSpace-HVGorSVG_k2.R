library(data.table)
library(Biostrings)
library(SpatialExperiment)
library(scater) # addPerCellQC
library(scran)
library(BayesSpace)
library(parallel)
library(BiocParallel)

#### Load list with SPEs of Harmony on top 10%ile, 20%ile HVGs and top 10%ile, 20%ile nnsvgs

pcas <- readRDS("analysis/data/spe_070123/H01-feature_selection/SPEs_sampleBlock-hvgs-10or20pctile_nnsvg10or20pctile_MNNk30-PCA-UMAP.RDS")

### note that reducedMNN returns the dimRed matrix as one column of a larger thing. so extract and append that to dimRed
pcas <- lapply(pcas,FUN=function(x){
    reducedDim(x,"mnn30.m") <- reducedDim(x,"mnn30")$corrected
    return(x)
})

### sanity check on the extracted reduced dims
stopifnot(sum(rownames(reducedDim(pcas[[1]],"mnn30.m"))==colnames(pcas[[1]]))==ncol(pcas[[1]]))

# kgo
bpparam <- MulticoreParam(4)
register(bpparam)
k2cls <- bplapply(pcas,BPPARAM=bpparam,FUN=function(x){
    ### according to https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/03_BayesSpace/01_BayesSpace.R , don't use spatial preprocess. in order to do this you have to reset metadata.
    ### HOWEVER, BayesSpace now throws an error if you don't actually run spatialpreproc. also, the documentation says that it ADDS metadata, which shouldn't be equiv to overwriting it. The "metadata" object in the HYP SPE is an empty list, so nothing will be modified. (Maybe LIBD used to place something in the metadata entry for something but evidently not now).
    colData(x)$row <- colData(x)$array_row
    colData(x)$col <- colData(x)$array_col
    x <- spatialPreprocess(x,platform="Visium",skip.PCA=T)


    # Main call
    # defaults to 50k reps, LIBD uses 10k
    # here, we're just using q=2 to separate out the two primary spatial components, whatever they may be.
    # tmp <- spatialCluster(sce = x,init.method = "mclust",use.dimred = "HARMONY",q = 2,platform = "Visium",nrep=2000)
    ## for spe070123, use mnn30
    tmp <- spatialCluster(sce = x,init.method = "mclust",use.dimred = "mnn30.m",q = 2,platform = "Visium",nrep=10000)


    # pull out clusters to save in much tinier file than an SPE
    tmp2 <- as.data.table(cbind(rownames(colData(tmp)),colData(tmp)$spatial.cluster))
    names(tmp2) <- c("rn","label")
    tmp2
})
names(k2cls) <- names(pcas)



i <- 1
for (i in c(1:4)){
write.table(k2cls[[i]],paste0("analysis/data/spe_070123/H02-clustering/03b-BSpace_mnn30_q2/",names(k2cls)[i],".txt"),sep='\t',quote=F,row.names=F,col.names=T)
rm(list=ls())
gc(full=T)
}

