# ml conda_R/4.3.x
# R

###########

dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20240223__174743__022324_KMon120823/"
# dir_data = "/dcs04/hansen/data/ywang/xenium/raw_Data/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_newRun/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out_QC_newRun/"
dir_data_processed_downstream = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"
dir.create(dir_data_processed_downstream)


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


####### MOLECULEEXPERIMENT to SpatialExperiment #######

######## load all Xenium results from their parent directory; change to understandable sample names; import cell (not nucleus) boundaries
xens <- MoleculeExperiment::readXenium(paste0(dir_data_processed,"/cellLevel/"),
                                       keepCols = "essential",
                                       addBoundaries = "cell")  


######## rename samples (which are the names of the molecules and boundaries parts of the molecule experiment)
# sample/section IDs by donor, spatial reps as letters; in order of slide #86 region 1-2-3, #97 reg 1-2-3 (how they will be loaded by default)
# xennames <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")
xennames <- c("br1225a","br1225b","br8741c",
              "br8741d","br5459a","br5459b",
              "br8667c")
sex = c("M", "M","F", 
        "F", "M", "M",
        "F")

for (i in c(1:length(xens@molecules$detected))){
  names(xens@molecules$detected)[i] <- xennames[i]
  names(xens@boundaries$cell)[i] <- xennames[i]
}



######## get names of features we want to keep (i.e. real genes)
keepgenes <- unlist(lapply(xens@molecules$detected,FUN=names))
keepgenes <- unique(keepgenes)
keepgenes <- grep(keepgenes,pattern="NegControlProbe|NegControlCodeword|DeprecatedCodeword|BLANK",value = T,invert = T)

######## test gene subsetting on two samples
######## subsetting MoleculeExperiment samples here doesn't have any straightforward approach off the top of my head (no function for it, anyhow), but i create a test object of just the first two samples to test the gene filtering. 
######## ends up being a useful **demonstration of how to create a new, identical moleculeexperimentobject with only samples of interest.**
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
# i<-1
# for (i in c(1:6)){
#   xens@molecules$detected[[i]] %<>% .[which(names(.) %in% keepgenes)]
# }
# rm(i)
# gc(full=T)

######## create a spatialexperiment for cell-level data ### 
######## with the molecularexperiment function **countMolecules**. specify cell boundaries as the assay we want to use. (This is infinitely more straightforward, and looks to be the most computationally efficient route I've identified, for building an SPE, including rather than trying to build a SPE manually/piecewise from Xenium output files).

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
saveRDS(xens, file=paste0(dir_data_processed_downstream, "xens.rds"))
