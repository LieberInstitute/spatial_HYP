library(data.table)
library(Biostrings)
library(SpatialExperiment)
library(scater) # addPerCellQC
library(scran)
library(BayesSpace)

#### Load list with SPEs of Harmony on top 10%ile, 20%ile HVGs and top 10%ile, 20%ile nnsvgs

pcas <- readRDS("analysis/data/H01-feature_selection/SPEs_sampleBlock-hvgs-10or20pctile_nnsvg10or20pctile_Harmony-PCA-UMAP.RDS")


### build table with constants (SPE to retrieve, by index in pcas, and value to pass to q in BayesSpace from among the optimized values of 9/15/20/31)



### pull params for task number


## important from slack: need to run the next two lines or will get "subscript contains invalid names"

bpparam <- MulticoreParam(4)
k2cls <- bplapply(pcas,BPPARAM=bpparam,FUN=function(x){
    ### according to https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/03_BayesSpace/01_BayesSpace.R , don't use spatial preprocess. in order to do this you have to reset metadata.
    ### HOWEVER, BayesSpace now throws an error if you don't actually run spatialpreproc. also, the documentation says that it ADDS metadata, which shouldn't be equiv to overwriting it. The "metadata" object in the HYP SPE is an empty list, so nothing will be modified. (Maybe LIBD used to place something in the metadata entry for something but evidently not now).
    colData(x)$row <- colData(x)$array_row
    colData(x)$col <- colData(x)$array_col
    spatialPreprocess(x,platform="Visium",skip.PCA=T)


    # Main call
    # defaults to 50k reps, LIBD uses 10k, we'll use 2k here just to see what separates first
    tmp <- spatialCluster(sce = x,init.method = "mclust",use.dimred = "HARMONY",q = 2,platform = "Visium",nrep=2000)

    # pull out clusters to save in much tinier file than an SPE
    tmp2 <- as.data.table(cbind(rownames(colData(tmp)),colData(tmp)$spatial.cluster))
    names(tmp2) <- c("rn","label")
    tmp2
})
names(k2cls) <- names(pcas)



i <- 1
for (i in c(1:4)){
write.table(k2cls[[i]],paste0("analysis/data/H02-clustering/03b-BSpace_Harmony_q2/",names(k2cls)[i],".txt"),sep='\t',quote=F,row.names=F,col.names=T)
rm(list=ls())
gc(full=T)
}

