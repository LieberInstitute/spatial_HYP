library(data.table)
library(Biostrings)
library(SpatialExperiment)
library(scater) # addPerCellQC
library(scran)
library(BayesSpace)

#### Load list with SPEs of Harmony on top 10%ile, 20%ile HVGs and top 10%ile, 20%ile nnsvgs

pcas <- readRDS("SPEs_sampleBlock-hvgs-10or20pctile_nnsvg10or20pctile_Harmony-PCA-UMAP.RDS")


### build table with constants (SPE to retrieve, by index in pcas, and value to pass to q in BayesSpace from among the optimized values of 9/15/20/31)

p <- rep(c(1:4),4)
q <- c(rep(9,4),rep(15,4),rep(20,4),rep(31,4))
pqs <- as.data.frame(cbind(p,q))
rm(p,q)


### pull params for task number

i <- as.numeric(Sys.getenv("SGE_TASK_ID"))
curdat <- pcas[[pqs[i,1]]]
curname <- names(pcas)[[pqs[i,1]]]
rm(pcas)


## important from slack: need to run the next two lines or will get "subscript contains invalid names"

colData(curdat)$row <- colData(curdat)$array_row
colData(curdat)$col <- colData(curdat)$array_col


### according to https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/03_BayesSpace/01_BayesSpace.R , don't use spatial preprocess. in order to do this you have to reset metadata.
### HOWEVER, BayesSpace now throws an error if you don't actually run spatialpreproc. also, the documentation says that it ADDS metadata, which shouldn't be equiv to overwriting it. The "metadata" object in the HYP SPE is an empty list, so nothing will be modified. (Maybe LIBD used to place something in the metadata entry for something but evidently not now).

spatialPreprocess(curdat,platform="Visium",skip.PCA=T)


# Main call
# defaults to 50k reps, LIBD uses 10k, we'll use 20k here

tmp <- spatialCluster(sce = curdat,init.method = "mclust",use.dimred = "HARMONY",q = pqs[i,2],platform = "Visium",nrep=20000)

# pull out clusters to save in much tinier file than an SPE
tmp2 <- as.data.table(cbind(rownames(colData(tmp)),colData(tmp)$spatial.cluster))
names(tmp2) <- c("rn",paste0("BShar_",curname,"_",pqs[i,2]))

write.table(tmp2,paste0("bs_harmony_9-15-20-31_out/BSpace_Harmony_",curname,"_",pqs[i,2],".txt"),sep='\t',quote=F,row.names=F,col.names=T)
rm(list=ls())
gc(full=T)
