# ml conda
# conda activate voyager
# R

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20240223__174743__022324_KMon120823/"
# dir_data = "/dcs04/hansen/data/ywang/xenium/raw_Data/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_newRun/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out_QC_newRun/"
dir_data_processed_downstream = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"
# dir.create(dir_data_processed_downstream)


###########
library(Voyager)
library(SFEData) 
library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(ggplot2)
library(stringr)
library(scater) 
library(scuttle)
library(BiocParallel)
library(BiocSingular)
library(bluster)
library(scran)
library(patchwork)
theme_set(theme_bw())

library(memisc)
library(here)
library(DropletUtils)
library(plotly)
library(arrow)
library(vroom)
library(nngeo)
library(listviewer)

######################################################################### 
######### load data
######################################################################### 
# load spe obj
# xens2 = readRDS(paste0(dir_data_processed, 
#                        "xens.rds"))
## function provided by Cindy:
readXenium <- function(dir_name){
  # Find the files
  counts_path <- here(dir_name, "cell_feature_matrix.h5")
  cell_info_path <- here(dir_name, "cells.csv.gz")
  cell_poly_path <- here(dir_name, "cell_boundaries.parquet")
  nuc_poly_path <- here(dir_name, "nucleus_boundaries.parquet")
  
  # Read in the data
  sce <- read10xCounts(counts_path)
  counts(sce) <- as(realize(counts(sce)), "dgCMatrix")
  cell_info <- vroom(cell_info_path)
  
  options(browser = 'firefox') # added by Yi
  
  cell_schema <- schema(cell_id=string(),
                        vertex_x=float64(),
                        vertex_y=float64())
  
  cell_poly_ <- open_dataset(cell_poly_path,
                             schema=cell_schema) 
  cell_poly = as.data.frame(cell_poly_)
  
  
  # cell_poly <- open_dataset(cell_poly_path,
  # schema=cell_schema) %>%
  # collect()
  
  nuc_poly <- open_dataset(nuc_poly_path,
                           schema=cell_schema) 
  # %>%
  # collect()
  nuc_poly = as.data.frame(nuc_poly)
  
  names(cell_poly)[1] <- "ID"
  names(nuc_poly)[1] <- "ID"
  
  cells_sf <- df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
  nuc_sf <- df2sf(nuc_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
  
  all(st_is_valid(cells_sf))
  all(st_is_valid(nuc_sf))
  
  ind_invalid <- !st_is_valid(nuc_sf)
  nuc_sf[ind_invalid,] <- nngeo::st_remove_holes(st_buffer(nuc_sf[ind_invalid,], 0))
  
  colData(sce) <- cbind(colData(sce), cell_info)
  print(cell_info)
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
  
  
  
  
  sfe <- toSpatialFeatureExperiment(spe)
  
  cellSeg(sfe, withDimnames = FALSE) <- cells_sf
  nucSeg(sfe, withDimnames = FALSE) <- nuc_sf
  
  
  # Add some QC Metrics
  colData(sfe)$nCounts <- colSums(counts(sfe))
  colData(sfe)$nGenes <- colSums(counts(sfe) > 0)
  
  is_blank <- str_detect(rownames(sfe), "^BLANK_")
  is_neg <- str_detect(rownames(sfe), "^NegControlProbe")
  is_neg2 <- str_detect(rownames(sfe), "^NegControlCodeword")
  is_anti <- str_detect(rownames(sfe), "^antisense")
  is_depr <- str_detect(rownames(sfe), "^DeprecatedCodeword")
  
  is_any_neg <- is_blank | is_neg | is_neg2 | is_anti | is_depr
  rowData(sfe)$is_neg <- is_any_neg
  
  n_panel <- nrow(sfe) - sum(is_any_neg)
  #print(n_panel)
  
  colData(sfe)$nCounts_normed <- sfe$nCounts/n_panel
  colData(sfe)$nGenes_normed <- sfe$nGenes/n_panel
  colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area
  
  sfe <- addPerCellQCMetrics(sfe, subsets = list(blank = is_blank,
                                                 negProbe = is_neg,
                                                 negCodeword = is_neg2,
                                                 anti = is_anti,
                                                 depr = is_depr,
                                                 any_neg = is_any_neg))
  
  
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  rowData(sfe)$cv2 <- rowData(sfe)$vars/rowData(sfe)$means^2
  
  
  # Add cell ids and make gene names unique
  colnames(sfe) <- seq_len(ncol(sfe))
  rownames(sfe) <- uniquifyFeatureNames(ID=rownames(sfe),  names=rowData(sfe)$Symbol)
  
  return(sfe)
}

######################################################################### 
# spe = xens2
dir_allSamples = dir_data
list_sample = list.files(dir_allSamples)
# length(list_sample)
# 6
for(sample_name in list_sample){
  print(sample_name)
  
  dir_name_ = paste0(dir_allSamples, sample_name)
  sfe_sample_tmp = readXenium(dir_name_)
  
  saveRDS(sfe_sample_tmp, file = paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
  
}

# dir_name = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20231019__161052__101923_KMon101923/output-XETG00089__00011586__Region_2__20231019__161214/"




######################################################################### 
# obj_sfe_oneSample = readXenium(dir_name)

# saveRDS(obj_sfe_oneSample, file = paste0(dir_data_processed, "obj_sfe_oneSample.rds"))





