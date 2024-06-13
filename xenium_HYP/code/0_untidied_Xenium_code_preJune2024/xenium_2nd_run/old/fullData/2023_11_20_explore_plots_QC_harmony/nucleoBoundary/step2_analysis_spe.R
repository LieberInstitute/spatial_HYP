# ml conda_R/4.3.x
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

######################################################################### 
######################################################################### 
##### start  from here!!!
######################################################################### 

xens = readRDS(paste0(dir_data_processed, "xens_nucleus.rds"))


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

xens2$brnum <- factor(xens2$brnum)
xens2$sample_id <- factor(xens2$sample_id)
xens2$slide <- factor(xens2$slide)


### gene variance
xenvars <- modelGeneVar(xens2,block=xens2$slide)
xens.vars.fit <- metadata(xenvars)
# get top 3/8ths of genes, since we have so few to work with
xenhvgs <- scran::getTopHVGs(xenvars,prop = 0.375)

### PCA's gonna need all genes though
library(Matrix)

xens2 <- scater::runPCA(xens2)
# xens2 <- scater::runPCA(xens2)
# Error in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  : 
                      # function 'as_cholmod_sparse' not provided by package 'Matrix'
# remotes::install_version("Matrix", version = "1.6-1") #
# it works after downgrading Matrix to 1.6-1 
saveRDS(xens2, file=paste0(dir_data_processed, 
                           "xens2_nucleus.rds"))#scater wasn;t run here

