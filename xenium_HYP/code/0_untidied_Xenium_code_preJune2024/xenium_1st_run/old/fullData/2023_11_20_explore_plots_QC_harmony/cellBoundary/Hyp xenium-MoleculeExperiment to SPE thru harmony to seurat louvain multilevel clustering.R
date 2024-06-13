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
library(MoleculeExperiment)
library(magrittr)
library(scater)
library(scran)
library(BiocParallel)
library(job)
library(SpatialExperiment)

###########
ColorOut()
options("styler.addins_style_transformer" = "biocthis::bioc_style()")

## change write.table to output a TSV without quotes or rownames by default
default(write.table) <- list(row.names=F,col.names=T,sep='\t',quote=F)

theme_set(theme_bw()+theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.y = element_text(size =16), plot.title = element_text(size = 20,hjust=0.5), strip.text = element_text(size=18), legend.text = element_text(size=10), legend.title = element_text(size=11,hjust=0.5)))


### Example seurat handling of the xenium data (up to the point of clustering, without any gene filters), is at the end of this Rmd. ###



####### MOLECULEEXPERIMENT to SpatialExperiment #######
# nucleus
## load all Xenium results from their parent directory; change to understandable sample names; import cell (not nucleus) boundaries
xens <- MoleculeExperiment::readXenium(paste0(dir_data,"/20231019__161052__101923_KMon101923/"),keepCols = "essential",addBoundaries = "cell")
# addBoundaries = "nucleus"


### rename samples (which are the names of the molecules and boundaries parts of the molecule experiment)
# sample/section IDs by donor, spatial reps as letters; in order of slide #86 region 1-2-3, #97 reg 1-2-3 (how they will be loaded by default)
xennames <- c("br6197","br5993a","br5993b",
              "br1735a","br1735b","br6588")
sex = c("M", "F", "F", "M", "M", "F")
# BrNum	AgeDeath	Race	Sex	PrimaryDx	PMI	Best RIN PFC	BMI (calculated)
# Br1735	31.54	CAUC	Male	Control	25	8.6	24.6
# Br5459	38.99	CAUC	Male	Control	32	8.5	23.1
# Br5993	48.15868	CAUC	Female	Control	28	6.2	27.1
# Br6197	28.55578	CAUC	Male	Control	42	9.1	24.5
# Br6588	27.34839	CAUC	Female	Control	31.5	6.2	27.3
# Br8406	30.95414	CAUC	Male	Control	28.5	6.8	22.3
# Br8667	37.33333	CAUC	Female	Control	13.5	6.9	31.6
# Br8741	28.19439	CAUC	Female	Control	21	7.4	24

i<-1
for (i in c(1:length(xens@molecules$detected))){
  names(xens@molecules$detected)[i] <- xennames[i]
  names(xens@boundaries$cell)[i] <- xennames[i]
}

### get names of features we want to keep (i.e. real genes)
keepgenes <- unlist(lapply(xens@molecules$detected,FUN=names))
keepgenes <- unique(keepgenes)
keepgenes <- grep(keepgenes,pattern="NegControlProbe|NegControlCodeword|DeprecatedCodeword|BLANK",value = T,invert = T)

### test gene subsetting on two samples
### subsetting MoleculeExperiment samples here doesn't have any straightforward approach off the top of my head (no function for it, anyhow), but i create a test object of just the first two samples to test the gene filtering. 
### ends up being a useful **demonstration of how to create a new, identical moleculeexperimentobject with only samples of interest.**
detectest <- list(xens@molecules$detected[[1]],xens@molecules$detected[[1]])
boundtest <- list(xens@boundaries$cell[[1]],xens@boundaries$cell[[2]])

names(detectest) <- xennames[c(1:2)]
names(boundtest) <- xennames[c(1:2)]
detectestholder <- list(detectest)
names(detectestholder) <- "detected"
boundtestholder <- list(boundtest)
names(boundtestholder) <- "cell"

test <- MoleculeExperiment(molecules=detectestholder,boundaries = boundtestholder)

i <-1
for (i in c(1:2)){
  test@molecules$detected[[i]] %<>% .[which(names(.) %in% keepgenes)]
}
stopifnot(length(names(test@molecules$detected[[1]]))==length(keepgenes))
## beautiful
rm(test,detectest,detectestholder,boundtest,boundtestholder,i)
gc(full=T)

### now really subset the moleculeexperiment
i<-1
for (i in c(1:6)){
  xens@molecules$detected[[i]] %<>% .[which(names(.) %in% keepgenes)]
}
rm(i)
gc(full=T)

### create a spatialexperiment for cell-level data ### 
### with the molecularexperiment function **countMolecules**. specify cell boundaries as the assay we want to use. (This is infinitely more straightforward, and looks to be the most computationally efficient route I've identified, for building an SPE, including rather than trying to build a SPE manually/piecewise from Xenium output files).

# This is able to perform multicore as long as BiocParallel is loaded.
xens <- countMolecules(xens,moleculesAssay = "detected",boundariesAssay = "cell",nCores = 8)

### holy SHIT that was fast. Like, 3 minutes.
gc(full=T)

# before we go further, we need to change the names of the "x_location" and "y_location" variables in the SPE's coldata, or ggspavis will get confused by that being names in both coldata and spatialcoords.
colnames(colData(xens))[3] <- "x_position"
colnames(colData(xens))[4] <- "y_position"

# save the resulting obj in the background
# job({saveRDS(xens,paste0(dir_data_processed, "hypxenia.cells.as.SPE_unfiltered.RDS"))#####start  from here!!!
  # job::export("none")},import="auto")
saveRDS(xens, file=paste0(dir_data_processed, "xens.rds"))

######################################################################### 
######################################################################### 
#####start  from here!!!
######################################################################### 




### now we can proceed with the usual approach
### gene-level mean-variance relationships


### get # of detections per cell and # of genes per cell, append to colData
detecs <- colSums(xens@assays@data$counts)
detecs <- DataFrame(detecs)
detecs <- detecs[colnames(xens),]
colData(xens) <- cbind(colData(xens),detecs)

ngene <- colSums(xens@assays@data$counts>0)
ngene <- DataFrame(as.data.frame(ngene))
ngene <- ngene[colnames(xens),]
colData(xens) <- cbind(colData(xens),ngene)
rm(ngene,detecs)
gc(full=T)
###



### filter to 20 counts / 10 genes

xens2 <- xens[,xens$detecs>19]
xens2 <- xens2[,xens2$ngene>9]
rm(xens)
gc(full=T)
## wow, only 1000 cells lost of 408k!

## calc norm factors, log normalize
xens2 <- scater::computeLibraryFactors(xens2)
xens2 <- scater::logNormCounts(xens2)


### add some additional variables
xens2$brnum <- gsub(xens2$sample_id, pattern="^(br[[:alnum:]]{4}).*$",replacement="\\1")
# make sure thsi worked unique(xens2$brnum)

# make a variable for slide
sldbrn <- as.data.table(as.data.frame(matrix(nrow=4,ncol=2)))
setnames(sldbrn,c("brnum","slide"))
sldbrn[,brnum:=unique(xens2$brnum)]
sldbrn[,slide:=as.character(slide)]
sldbrn[brnum %in% c("br6197","br5993"),slide:="slide86"]
sldbrn[brnum %in% c("br1735","br6588"),slide:="slide97"]

tmpcd <- as.data.table(colData(xens2),keep.rownames=T)
tmpcd <- merge.data.table(tmpcd,sldbrn,by="brnum")
tmpcd <- DataFrame(tmpcd)
rownames(tmpcd) <- tmpcd$rn
tmpcd$rn <- NULL
tmpcd <- tmpcd[colnames(xens2),]
colData(xens2) <- tmpcd
rm(tmpcd,sldbrn)

## make sure the variables (sample, donor, slide) are factors
xens2$brnum <- factor(xens2$brnum)
xens2$sample_id <- factor(xens2$sample_id)
xens2$slide <- factor(xens2$slide)


### gene variance
xenvars <- modelGeneVar(xens2,block=xens2$slide)
xens.vars.fit <- metadata(xenvars)
# get top 3/8ths of genes, since we have so few to work with
xenhvgs <- scran::getTopHVGs(xenvars,prop = 0.375)

### PCA's gonna need all genes though
xens2 <- scater::runPCA(xens2)

# sbp <- MulticoreParam(6)
# cluss <- quickCluster(x = xens2,block=xens2$sample_id,BPPARAM=sbp)
# ### lol what? 125 clusters?
# clus2 <- cbind(colnames(xens2),cluss)
# clus2 <- DataFrame(clus2)
# 
# saveRDS(clus2, "~/Desktop/quickcluss_sampleblock_xens.RDS")
# 
# colLabels(xens2) <- as.factor(paste0("X",clus2$cluss))
rm(clus2,xenvars,xens.vars.fit)
gc(full=T)
# plot that shit!

# clus125 <- Polychrome::createPalette(125,seedcolors = Polychrome::palette36.colors(),range = c(30,90))
# names(clus125) <- paste0("X",c(1:125))
# 
# plt <- ggspavis::plotSpots(xens2,sample_id = "sample_id",annotate = "label",size = 0.2,palette = clus125,in_tissue = NULL)
# 
# dev.off()
# pdf("~/Desktop/alottacells.pdf",height=100,width=100)
# plt
# dev.off()



### ^ okay, but 125 clusters doesn't even make remote sense, so
library(BiocParallel)
sbp <- MulticoreParam(2) # this will max out at 2 processes spawned even if we specify more here, since we're blocking on a factor with 2 levels. and even just spawning two processes here is going to need 100 GB of RAM+swap for k=10

# k=10 walktrap
cluss <- quickCluster(x = xens2,block=xens2$slide,BPPARAM=sbp,method="igraph",graph.fun="walktrap",k=10)
### ?
clus2 <- cbind(colnames(xens2),cluss)
clus2 <- DataFrame(clus2)
rownames(clus2) <- clus2[,1]
colnames(clus2) <- c("cellid","cluster")
saveRDS(clus2,"~/Desktop/xens_slideBlock_walktrap-k10_quickclus.RDS")

clus2$cluster <- paste0("X",clus2$cluster)
clus2 <- clus2[colnames(xens2),]
colLabels(xens2) <- factor(clus2$cluster)

gc(full=T)
# 
# clus51 <- Polychrome::createPalette(51,seedcolors = Polychrome::palette36.colors()[3:36],range = c(30,90))
# names(clus51) <- paste0("X",c(1:51))
# # 
# plt <- ggspavis::plotSpots(xens2,sample_id = "sample_id",annotate = "label",size = 0.2,palette = clus51,in_tissue = NULL)
# 
# pdf("~/Desktop/alottacellsinfewerclusters.pdf",height=100,width=100)
# plt
# dev.off()

### still assigning different clusters by sex. that won't do. prob because XIST (definitely sex DE) is in there driving variation.


### HARMONY
harmo <- harmony::RunHarmony(xens2,group.by.vars=c("brnum","slide"),ncores=3)
rm(xens2)
gc(full=T)
## save results
job::job({saveRDS(harmo,"~/Desktop/hypxenia.cells.as.SPE_filtered_Harmony-brnum-slide.RDS")
  job::export("none")},import = "auto")
##


## seurat louvain multilevel refinement clustering? 
xens3 <- as.Seurat(xens2)

library(future)
plan(multisession,workers=6)
plan()

cluss <- Seurat::FindNeighbors(xens3,reduction="HARMONY",dims=1:30,k.param = 30)

cluss <- Seurat::FindClusters(cluss,method="igraph",algorithm = 2,random.seed = 42,verbose = T)

### THERE we go. 35 clusters. (note that algo 2, louvain meltivel, was the most consistently strong across with semi-synthetic spatial xscripto datasets in a review of clustering methods, even though its spatially unaware). incidentally, seurat's approaches were all the best performers for the one merfish dataset included in that study.

seucls <- Idents(cluss)
seucls <- DataFrame(as.data.frame(seucls))
seucls <- as.data.table(seucls,keep.rownames=T)
setnames(seucls,2,"SPEHarmony_seurat_louvMultiRefine_clus")
### setting these as numeric for the moment also, conveniently, brought them to numbering that counts from 1 rather than 0.
seucls[,SPEHarmony_seurat_louvMultiRefine_clus:=as.numeric(SPEHarmony_seurat_louvMultiRefine_clus)]
# fwrite(seucls,"~/Desktop/23-10-24-Xen/xenium sandbox/data/HYPxen.filt.Harmony-brnum-slide.seuratK30SNN-louvainMultilevelRefine_clusters.txt",sep='\t',quote=F)
rm(xens3,cluss)
unloadNamespace("Seurat")
gc(full=T)



### examine by plotting: post-harmony louvain clusters
seucls <- DataFrame(seucls)
rownames(seucls) <- seucls$rn
seucls$SPEHarmony_seurat_louvMultiRefine_clus <- paste0("X",seucls$SPEHarmony_seurat_louvMultiRefine_clus)
seucls <- seucls[colnames(xens2),]

colLabels(xens2) <- as.factor(seucls$SPEHarmony_seurat_louvMultiRefine_clus)


library(ggspavis)
clus35 <- scCustomize::DiscretePalette_scCustomize(num_colors = 35,palette = "varibow",shuffle_pal = F) # this won't have any blacks  or grays, which are two colors of the polychrome 36-color, annoyingly

names(clus35) <- paste0("X",c(1:35))
plt <- ggspavis::plotSpots(xens2,annotate="label",palette = clus35,in_tissue = NULL,size = 0.3)

pdf("~/Desktop/harmonybrnumslide_seuratlouvaink30clus.pdf",height=100,width=100)
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
anchors <- FindIntegrationAnchors(xenlist,anchor.features = features,reference=c(1,3,4,6),reduction="rpca",dims=1:50,normalization.method = "SCT",k.filter = 150,verbose = T)

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




