# ml conda_R/4.3
# R


###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out/"

###########
library(data.table)
library(Biostrings)
library(ggplot2)
library(gridExtra)
#library(SpatialExperiment)
library(spatialLIBD)
require(colorout)
library(default)
# library(MoleculeExperiment)
library(magrittr)
library(scater)
library(scran)
library(BiocParallel)
library(job)
library(SpatialExperiment)
library(BiocParallel)
library(ggspavis)
library(Seurat)

###########
ColorOut()
options("styler.addins_style_transformer" = "biocthis::bioc_style()")

## change write.table to output a TSV without quotes or rownames by default
default(write.table) <- list(row.names=F,col.names=T,sep='\t',quote=F)

theme_set(theme_bw()+theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.y = element_text(size =16), plot.title = element_text(size = 20,hjust=0.5), strip.text = element_text(size=18), legend.text = element_text(size=10), legend.title = element_text(size=11,hjust=0.5)))

######################################################################### 

xens2 = readRDS(paste0(dir_data_processed, 
                       "xens2.rds"))
######################################################################### 
######################################################################### 

### HARMONY
head(colData(xens2))
table(colData(xens2)$brnum)
table(colData(xens2)$slide)

harmo <- harmony::RunHarmony(xens2,lambda=NULL,group.by.vars=c("brnum","slide"),ncores=1)
# rm(xens2)
# gc(full=T)

##
saveRDS(harmo, file = paste0(dir_data_processed, 
                             "xens2_harmony.rds"))
harmo = readRDS( paste0(dir_data_processed, 
                        "xens2_harmony.rds"))


## seurat louvain multilevel refinement clustering? 
# xens3 <- as.Seurat(xens2) # deleted by Yi
xens3  <- as.Seurat(harmo) # added by Yi
# saveRDS(xens3, file = paste0(dir_data_processed, 
                             # "xens3.rds"))
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from x_location to xlocation_
# Warning: All keys should be one or more alphanumeric characters followed by an underscore '_', setting key to xlocation_
# Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from PC to PC_
# Warning: All keys should be one or more alphanumeric characters followed by an underscore '_', setting key to PC_


cluss <- Seurat::FindNeighbors(xens3,
                               reduction="HARMONY",
                               dims=1:30,k.param = 30)
saveRDS(cluss, file = paste0(dir_data_processed,
  "cluss_NN.rds"))

cluss <- Seurat::FindClusters(cluss,method="igraph",algorithm = 2,random.seed = 42,verbose = T)
saveRDS(cluss, file = paste0(dir_data_processed,
                             "cluss.rds"))


### THERE we go. 35 clusters. (note that algo 2, louvain meltivel, was the most consistently strong across with semi-synthetic spatial xscripto datasets in a review of clustering methods, even though its spatially unaware). incidentally, seurat's approaches were all the best performers for the one merfish dataset included in that study.
seucls <- Idents(cluss)
saveRDS(seucls, file = paste0(dir_data_processed,
                             "idents_seucls.rds"))

##############################################################
### plot - start from here
##############################################################
xens2 = readRDS(paste0(dir_data_processed, 
                       "xens2.rds"))
# cluss = readRDS(paste0(dir_data_processed,
                       # "cluss.rds"))
cluss = readRDS(paste0(dir_data_processed,
                        "cluss.rds"))
seucls = readRDS(paste0(dir_data_processed,
                        "idents_seucls.rds"))
table(seucls) # 37 clusters
seucls <- DataFrame(as.data.frame(seucls))
seucls <- as.data.table(seucls, keep.rownames=T)
setnames(seucls,2,"SPEHarmony_seurat_louvMultiRefine_clus")
### setting these as numeric for the moment also, conveniently, brought them to numbering that counts from 1 rather than 0.
seucls[,SPEHarmony_seurat_louvMultiRefine_clus:=as.numeric(SPEHarmony_seurat_louvMultiRefine_clus)]
# fwrite(seucls,"~/Desktop/23-10-24-Xen/xenium sandbox/data/HYPxen.filt.Harmony-brnum-slide.seuratK30SNN-louvainMultilevelRefine_clusters.txt",sep='\t',quote=F)
# rm(xens3,cluss)
# unloadNamespace("Seurat")
# gc(full=T)




### examine by plotting: post-harmony louvain clusters
seucls <- DataFrame(seucls)
rownames(seucls) <- seucls$rn
seucls$SPEHarmony_seurat_louvMultiRefine_clus <- paste0("X",seucls$SPEHarmony_seurat_louvMultiRefine_clus)
seucls <- seucls[colnames(xens2),]

colLabels(xens2) <- as.factor(seucls$SPEHarmony_seurat_louvMultiRefine_clus)




clus35 <- scCustomize::DiscretePalette_scCustomize(num_colors = 37,palette = "varibow",shuffle_pal = F) # this won't have any blacks  or grays, which are two colors of the polychrome 36-color, annoyingly
# install.packages("scCustomize")
# Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) : 
#   namespace ‘SeuratObject’ 4.1.4 is already loaded, but >= 5.0.0 is required
# Calls: <Anonymous> ... namespaceImportFrom -> asNamespace -> loadNamespace
# Execution halted
# ERROR: lazy loading failed for package ‘scCustomize’
# * removing ‘/users/ywang/R/4.3/scCustomize’
# 
# The downloaded source packages are in
# ‘/tmp/Rtmp6fJhff/downloaded_packages’
# Warning messages:
#   1: In install.packages("scCustomize") :
#   installation of package ‘SeuratObject’ had non-zero exit status
# 2: In install.packages("scCustomize") :
#   installation of package ‘scCustomize’ had non-zero exit status

## install older version - worked
# install.packages(scCustomize_1.1.0.tar.gz")

names(clus35) <- paste0("X",c(1:37))
saveRDS(clus35, file = paste0(dir_out,"pallete_37clus.rds"))

plt <- ggspavis::plotSpots(xens2, # annotated by Yi
                           annotate="label", 
                           # palette="label", #added by Yi
                           palette = clus35, # annotated by Yi
                           in_tissue = NULL,
                           size = 0.3)

pdf(paste0(dir_out, "harmonybrnumslide_seuratlouvaink30clus.pdf"),height=100,width=100)
plt
dev.off()






#####################
#  Xenium Handling  #
####  w Seurat:  ####
#####  Example  #####
#####################



### EXAMPLE--NOT USED### 
### SEURATv5 Xenium pipeline: single sample example ###
# library(Seurat)
# xenium.obj <- Seurat::loadXenium("/Users/bmulvey/Desktop/23-10-24-Xen/20231019__161052__101923_KMon101923/output-XETG00089__00011597__Region_1__20231019__161215")
# xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
# xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
# 
# ### SCTransform actually defines variable features already, so they don't need to be specified for PCA. (We have to select integration features and use those only in the case of multi-sample)
# xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
# xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
# xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
# xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)


### EXAMPLE--NOT USED###
### Seurat xenium, all samples - SCtransform normalization - rpca integration ###  
library(BiocParallel)
sbp <- MulticoreParam(6)
register(sbp)

xens <- paste0(list.files ("/Users/bmulvey/Desktop/23-10-24-Xen/20231019__161052__101923_KMon101923",full.names = T),"/")
xenlist <- bplapply(xens,FUN=LoadXenium,BPPARAM = sbp)
rm(sbp)
gc(full=T)

#### seurat can parallelize with future::, so unload biocparallel in the meantime to avoid conflicts
unloadNamespace("BiocParallel")
####

# sample/section IDs by donor, spatial reps as letters; in order of slide #86 region 1-2-3, #97 reg 1-2-3 (how they will be loaded by default)
names(xenlist) <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

### get rid of control count assays (which are stored separately from the gene counts by seurat)
i<-1
for (i in c(1:6)){
  xenlist[[i]]@assays$BlankCodeword <- NULL
  xenlist[[i]]@assays$ControlCodeword <- NULL
  xenlist[[i]]@assays$ControlProbe <- NULL
  xenlist[[i]]$donor_id <- names(xenlist)[i]
}
rm(i)
gc(full=T)

###### set up future for parallelized computing
library(future)
library(future.apply)
plan(multisession,workers=6)
plan()
###### set max RAM that can be allocated to a parallelized process AT ITS START (needs to be > size of packages+data object(s)) 
options(future.globals.maxSize = 11010048000)


#### Normalize/scale/variable gene detection using Seurat's SCTransform
xenlist <- future_lapply(X = xenlist, future.seed=TRUE,FUN = function(X){ 
  Y <- subset(X,subset=nCount_Xenium>0)
  Y <- SCTransform(Y,assay="Xenium",clip.range=c(-10,10),variable.features.n=274)#top 75% of variable features per sample --> consensus 37.5% will be determined across all samples at next step
  return(Y)}
)
#### make sure parallel sessions stop (specifically to free held memory)
plan(sequential)
plan(multisession,workers=6)

#### choose top 37.5% of genes as genes for integration across samples
features <- SelectIntegrationFeatures(xenlist,nfeatures=137)

plan(sequential)
plan(multisession,workers=6)

#### perform initial PCA using the integration genes
xenlist <- future_lapply(X=xenlist,future.seed=TRUE,FUN=function(X){
  Y <- RunPCA(X,features=features)
  return(Y)}
)

plan(sequential)
gc(full=T)
plan(multisession,workers=6)

#### necessary but unsure what this step does
xenlist <- PrepSCTIntegration(xenlist,anchor.features = features,verbose = T)

plan(sequential)
gc(full=T)
plan(multisession,workers=4)

# use each sex per slide (1 M, 1 F)x2 as the "reference" samples for integration -- i.e., spannign variability across slides and sexes, ideally. the references are integrated in the order (a,b,c,d) -->
# 1. C+D
# 2. B + (CD)
# 3. A + (BCD)
# so we want to make sure this list is something like (slide1male,slide2female,slide1male,slide2male) so that we integrate accounting for the technical, rather than biological, variable of interest first. im not quite positive about that so DO YOUR OWN READING BEFORE TAKING THIS FOR FACT.

### regardless, we don't reach the finish line here because we'll run out of memory on a 32 GB RAM+600GB swap when trying to actually integrate a couple calls down from here.
anchors <- FindIntegrationAnchors(xenlist,
                                  anchor.features = features,
                                  reference=c(1,3,4,6),
                                  lambda=NULL,
                                  reduction="rpca",
                                  dims=1:50,
                                  normalization.method = "SCT",
                                  k.filter = 150,
                                  verbose = T)

plan(sequential)
gc(full=T)

### we will use seurat's "rpca" (reciprocal PCA) sample integration approach, which is supposedly streamlined for large datasets. supposedly.
plan(multisession,workers=2)
options(future.globals.maxSize = 26010048000) 
# ^ increase to approx 25 GB for the integration step - the partially integrated objects are returned to the environment and then sent back out to future processes after each sample integration (which grows the integrated data in progress, and means the outgoing partially-integrated data has to fall under globals.maxsize too.)


### NVM - this next step runs out of memory (i.e. depletes, swap, more than 500GB) before completing.

### NOT RUN (successfully): ###
# xeninteg <- IntegrateData(anchorset = anchors,dims=1:50,normalization.method = "SCT",verbose = T)

# saveRDS(xeninteg,file="~/Desktop/xen_hyp_integ_SCT-rpca.RDS")




