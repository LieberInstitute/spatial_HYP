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
