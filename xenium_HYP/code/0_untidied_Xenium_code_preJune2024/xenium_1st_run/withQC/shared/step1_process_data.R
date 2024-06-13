# ml conda_R/4.3.x
# R

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20231019__161052__101923_KMon101923/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out_QC_nuclear/"
dir_data_processed_downstream = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_nuclear/"
# dir.create(dir_data_processed_downstream)
# dir.create(dir_out)

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


####### MOLECULEEXPERIMENT to SpatialExperiment #######

######## load all Xenium results from their parent directory; change to understandable sample names; import cell (not nucleus) boundaries
xens <- MoleculeExperiment::readXenium(paste0(dir_data_processed,"/nucleusLevel/"),
                                       keepCols = "essential",
                                       addBoundaries = "nucleus")  


######## rename samples (which are the names of the molecules and boundaries parts of the molecule experiment)
# sample/section IDs by donor, spatial reps as letters; in order of slide #86 region 1-2-3, #97 reg 1-2-3 (how they will be loaded by default)
xennames <- c("br1225a","br1225b","br8741c",
              "br8741d","br5459a","br5459b",
              "br8667c")

for (i in c(1:length(xens@molecules$detected))){
  names(xens@molecules$detected)[i] <- xennames[i]
  names(xens@boundaries$nucleus)[i] <- xennames[i]
}

ls(names(xens@molecules$detected)[[1]])
# > names(xens@molecules$detected[[1]])
# [1] "ABCC9"                   "ACVR1C"                 
# [3] "ADAMTS12"                "ADAMTS16"               
# [29] "BDNF"                    "BLANK_0001"             
# [31] "BLANK_0005"              "BLANK_0007"             
# [103] "BLANK_0213"              "BLANK_0217"             
# [143] "CALCRL"                  "CAMK2A"                 
# [145] "CAPG"                    "CAPN3"                  
# [203] "DLG4"                    "DNER"                   
# [205] "DOC2A"                   "DYDC2"                  
# [207] "DeprecatedCodeword_0160" "DeprecatedCodeword_0171"
# [213] "DeprecatedCodeword_0373" "ECEL1"                  
# [215] "EFHD1"                   "EGFR"                   
# [341] "NRP1"                    "NTNG1"                  
# [345] "NWD2"                    "NXPH2"                  
# [399] "NegControlProbe_00024"   "NegControlProbe_00025"  
# [405] "NegControlProbe_00039"   "NegControlProbe_00041"  
# [407] "NegControlProbe_00042"   "OLIG1"                  
# [409] "OLIG2"                   "OPALIN"                 
# [411] "OTOGL"                   "OXTR"                   
# [423] "PCSK6"                   "PDGFD"                  
# [425] "PDGFRA"                  "PECAM1"                 


# ######## get names of features we want to keep (i.e. real genes)
# keepgenes <- unlist(lapply(xens@molecules$detected,FUN=names))
# keepgenes <- unique(keepgenes)
# keepgenes <- grep(keepgenes,pattern="NegControlProbe|NegControlCodeword|DeprecatedCodeword|BLANK",value = T,invert = T)
# 
# ######## test gene subsetting on two samples
# ######## subsetting MoleculeExperiment samples here doesn't have any straightforward approach off the top of my head (no function for it, anyhow), but i create a test object of just the first two samples to test the gene filtering. 
# ######## ends up being a useful **demonstration of how to create a new, identical moleculeexperimentobject with only samples of interest.**
# detectest <- list(xens@molecules$detected[[1]],xens@molecules$detected[[1]])
# boundtest <- list(xens@boundaries$cell[[1]],xens@boundaries$cell[[2]])
# 
# names(detectest) <- xennames[c(1:2)]
# names(boundtest) <- xennames[c(1:2)]
# detectestholder <- list(detectest)
# names(detectestholder) <- "detected"
# boundtestholder <- list(boundtest)
# names(boundtestholder) <- "cell"
# 
# test <- MoleculeExperiment(molecules=detectestholder,boundaries = boundtestholder)
# 
# i <-1
# for (i in c(1:2)){
#   test@molecules$detected[[i]] %<>% .[which(names(.) %in% keepgenes)]
# }
# stopifnot(length(names(test@molecules$detected[[1]]))==length(keepgenes))
## beautiful
# rm(test,detectest,detectestholder,boundtest,boundtestholder,i)
# gc(full=T)
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
xens <- countMolecules(xens,moleculesAssay = "detected",boundariesAssay = "nucleus",nCores = 8)

### holy SHIT that was fast. Like, 3 minutes.
gc(full=T)

# before we go further, we need to change the names of the "x_location" and "y_location" variables in the SPE's coldata, or ggspavis will get confused by that being names in both coldata and spatialcoords.
colnames(colData(xens))[3] <- "x_position"
colnames(colData(xens))[4] <- "y_position"

# save the resulting obj in the background
# job({saveRDS(xens,paste0(dir_data_processed, "hypxenia.cells.as.SPE_unfiltered.RDS"))#####start  from here!!!
  # job::export("none")},import="auto")
saveRDS(xens, file=paste0(dir_data_processed_downstream, "xens.rds"))
