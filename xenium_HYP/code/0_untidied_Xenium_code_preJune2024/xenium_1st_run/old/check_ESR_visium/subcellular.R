# ml conda_R/4.3.x
# R
# https://sydneybiox.github.io/MoleculeExperiment/articles/MoleculeExperiment.html
library(MoleculeExperiment)
library(ggplot2)

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
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out/"
###########
xens <- MoleculeExperiment::readXenium(paste0(dir_data,"/20231019__161052__101923_KMon101923/"),keepCols = "essential",addBoundaries = "cell")

tmpxens2 <- (xens)
tmpxens2@molecules <- lapply(tmpxens2@molecules,FUN=function(x){
  y <- lapply(x,FUN=function(x){
    z <- c("ESR1","ESR1") # or your vector of gene names of interest here
    return(z)
  })
  return(y)
})

### sub-subset to minimum of 2 samples (can replace indices with sample ids eg “br6197”)
tmpxens2@molecules$detected <- tmpxens2@molecules$detected
tmpxens2@boundaries$cell <- tmpxens2@boundaries$cell

ls(tmpxens2@molecules$detected[[1]])

### pass to moleculeexperiment plotting functions
# ggplot_me() +
#   #geom_polygon_me(xens, byFill = “cell”, colour = “black”) +
#   # geom_point_me(tmpxens2, byColour = "feature_id", size = 0.1)+
#   geom_point_me(tmpxens2,size = 0.1)+
#   geom_polygon_me(tmpxens2, assayName = "cell", fill = NA, colour = "red")

pdf("esr1_subcellular.pdf")
ggplot_me() +
  geom_polygon_me(tmpxens2, assayName = "cell", fill = "grey") +
  geom_point_me(tmpxens2) #+
  # zoom in to selected patch area
  # coord_cartesian(
  #   xlim = c(4900, 4919.98),
  #   ylim = c(6400.02, 6420)
  # )

dev.off()




