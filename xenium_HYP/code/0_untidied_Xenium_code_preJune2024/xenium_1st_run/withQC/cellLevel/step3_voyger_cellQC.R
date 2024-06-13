# ml conda
# conda activate voyager
# R

# based on this tutorial: https://pachterlab.github.io/voyager/articles/vig5_xenium.html

###########
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20231019__161052__101923_KMon101923/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out_QC/"
dir_data_processed_downstream = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
setwd(dir_out)
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

theme_set(theme_bw())


library(cowplot)
# library(RColorBrewer)
library(grid)
######################################################################### 
######### load data
######################################################################### 
# load spe obj
# xens2 = readRDS(paste0(dir_data_processed_downstream, 
# "xens2.rds"))

######################################################################### 
dir_allSamples = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20231019__161052__101923_KMon101923/"
list_sample = list.files(dir_allSamples)

######################################################################### 


# sample_name = list_sample[1]
# sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))

plist = list()
pdf("cellSeg.pdf")
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
  colnames(sfe) <- seq_len(ncol(sfe))
  plist[[ksample]] = plotGeometry(sfe, "cellSeg")
}
print(plot_grid(plotlist = plist,nrow=2))
dev.off()



make_plot_grid<- function(fun_voyager, name_pdf, main, list_arg, 
                          if_sfe = TRUE, 
                          name_arg_first = NULL,
                          sfe_keep = FALSE){
  # if_sfe: if the first argument of fun_voyager is 'sfe'
  plist = list()
  pdf(name_pdf, width = 30, height = 20)
  for(ksample in 1:length(list_sample)){
    sample_name = list_sample[ksample]
    if(sfe_keep){
      sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_keep_oneSample_",sample_name,".rds"))
    }else{
      sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
    }
    colnames(sfe) <- seq_len(ncol(sfe))
    
    if(if_sfe){
      list_arg$sfe = sfe
    }else{
      list_arg[[name_arg_first]] = sfe
    }
    plist[[ksample]] = do.call(fun_voyager, list_arg)  + ggtitle(paste0(main,", sample ", ksample))
  }
  print(plot_grid(plotlist = plist,nrow=2))
  dev.off()
  
}


make_plot_grid_plist <- function(fun_voyager, main, list_arg, 
                                 if_sfe = TRUE,
                                 name_arg_first = NULL, 
                                 sfe_keep= FALSE){
  # if_sfe: if the first argument of fun_voyager is 'sfe'
  plist = list()
  # pdf(name_pdf, width = 30, height = 20)
  for(ksample in 1:length(list_sample)){
    sample_name = list_sample[ksample]
    if(sfe_keep){
      sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_keep_oneSample_",sample_name,".rds"))
    }else{
      sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
    }
    colnames(sfe) <- seq_len(ncol(sfe))
    
    if(if_sfe){
      list_arg$sfe = sfe
    }else{
      list_arg[[name_arg_first]] = sfe
    }
    plist[[ksample]] = do.call(fun_voyager, list_arg)  + ggtitle(paste0(main,", sample ", ksample))
  }
  # print(plot_grid(plotlist = plist,nrow=2))
  # dev.off()
  
  plist
}

# This is what the tissue, with the cell outlines, looks like
make_plot_grid(fun_voyager = plotGeometry, 
               name_pdf = "cellSeg.pdf",
               main = "cell segmentation",
               list_arg = list(type="cellSeg"))

# Plot cell density in space
make_plot_grid(fun_voyager = plotCellBin2D, 
               name_pdf = "cellDensityInSpace.pdf",
               main = "cell density in space",
               list_arg = list(hex=TRUE))

###### Quality control 
#### Cells
# Some QC metrics are precomputed and are stored in colData
# names(colData(sfe))


# Since there’re more cells, it would be better to plot the tissue larger, so we’ll plot the histogram 
# of QC metrics and the spatial plots separately, unlike in the CosMx vignette.

n_panel <- 313
colData(sfe)$nCounts_normed <- sfe$nCounts/n_panel
colData(sfe)$nGenes_normed <- sfe$nGenes/n_panel

# Here we divided nCounts by the total number of genes probed, so this histogram is comparable to those from other smFISH-based datasets.
make_plot_grid(fun_voyager = plotColDataHistogram, 
               name_pdf = "density_nCounts_nGene_histogram.pdf",
               main = "nCounts, nGene",
               list_arg = list(feature=c("nCounts_normed", "nGenes_normed")))

 
# Compared to the FFPE CosMX non-small cell lung cancer dataset, 
# more transcripts per gene on average and a larger proportion of all genes are 
# detected in this dataset, which is also FFPE. However, this should be interpreted with care, 
# since these two datasets are from different tissues and have different gene panels, 
# so this may or may not indicate that Xenium has better detection efficiency than CosMX.
make_plot_grid(fun_voyager = plotSpatialFeature, 
               name_pdf = "nCounts_inSpace.pdf",
               main = "nCounts in space",
               list_arg = list(feature=c("nCounts"), 
                               colGeometryName = "cellSeg",
                               color = NA))# the cell boundaries makes each dot looks so dark, hard to see the blue color inside

## additional plot that Stephanie is interested
make_plot_grid(fun_voyager = plotColData, 
               name_pdf = "nCounts_cellArea_scatter.pdf",
               main = "cell area vs nCounts",
               list_arg = list(x = "cell_area",
                               y = "nCounts"),
               if_sfe = FALSE,
               name_arg_first = "object")

# color = "black"
# There seem to be FOV artifacts. However, the cell ID
# and FOV information were unavailable so we cannot examine them.
make_plot_grid(fun_voyager = plotSpatialFeature, 
               name_pdf = "nGenes_inSpace.pdf",
               main = "nGenes in space",
               list_arg = list(feature=c("nGenes"), 
                               colGeometryName = "cellSeg",
                               color= NA)) # the cell boundaries makes each dot looks so dark, hard to see the blue color inside

# A standard examination is to look at the relationship between nCounts and nGenes:
make_plot_grid(fun_voyager = plotColData, 
               name_pdf = "nGenes_vs_nCounts_scatter.pdf",
               main = "nGenes vs nCounts",
               list_arg = list(x = "nCounts", 
                               y = "nGenes",
                               bins = 100),
               if_sfe = FALSE,
               name_arg_first = "object")



# There appear to be two branches.
# Here we plot the distribution of cell area
# plotColDataHistogram(sfe, c("cell_area", "nucleus_area"), scales = "free_y")
make_plot_grid(fun_voyager = plotColDataHistogram, 
               name_pdf = "density_cellArea_nucleusArea.pdf",
               main = "cell area, nucleus area",
               list_arg = list(feature = c("cell_area", "nucleus_area"),
                               scales = "free_y"))

# That should be in pixels. There’s a very long tail. 
# The nuclei are much smaller than the cells.
# How is cell area distributed in space?
make_plot_grid(fun_voyager = plotSpatialFeature, 
               name_pdf = "cellArea_inSpace.pdf",
               main = "cell area distributed in space",
               list_arg = list(feature = c("cell_area"),
                               scales = "cellSeg",
                               color = NA))


# Cells in the sparse region tend to be larger than those in the dense region.
# This may be biological or an artifact of the cell segmentation algorithm or both.

# Here the nuclei segmentations are plotted instead of cell segmentation. 
# The nuclei are much smaller to the extent that they are difficult to see.
make_plot_grid(fun_voyager = plotSpatialFeature, 
               name_pdf = "nucleusArea_inSpace.pdf",
               main = "nucleus area distributed in space",
               list_arg = list(feature = c("nucleus_area"),
                               colGeometryName = "nucSeg",
                               color = NA))


# There’s an outlier near the right edge of the section, throwing off 
# the dynamic range of the plot. Upon inspection of the H&E image, 
# the outlier is a bit of tissue debris that doesn’t look like a cell. 
# But we can still that cells in the dense, gland like regions tend to 
# have larger nuclei. This may be biological, or that nuclei are so densely 
# packed in those regions that they are more likely to be undersegmented, 
# i.e. when multiple nuclei are counted as one by the nuclei segmentation 
# program, or both.

# These observations motivate an examination of the relationship between 
# cell area and nuclei area:

# plotColData(sfe, x="cell_area", y="nucleus_area", bins = 100)

make_plot_grid(fun_voyager = plotColData, 
               name_pdf = "nucleusArea_cellArea_scatter.pdf",
               main = "cell area vs. nucleus area",
               list_arg = list(x = "cell_area",
                               y = "nucleus_area", 
                               bins = 100),
               if_sfe = FALSE,
               name_arg_first = "object")


# Again, there are two branches, probably related to cell density and cell type. 
# The nucleus outlier also has large cell area, though it is not as much an outlier 
# in cell area. However, it is a spatial outlier as it’s unusually large compared to 
# its neighbors (scroll up two plots back).

# Next we calculate the proportion of cell in this z-plane taken up by the nucleus, 
# and examine the distribution:
# Add prop_nuc to the colData of sfe file and save, for each sample
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
  # colnames(sfe) <- seq_len(ncol(sfe))
  colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area
  # save the updated sfe file with new colData column added
  saveRDS(sfe, file = paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
}
make_plot_grid(fun_voyager = plotColDataHistogram, 
               name_pdf = "prop_nuc_histgram.pdf",
               main = " proportion of cell in this z-plane taken up by the nucleus",
               list_arg = list(feature = "prop_nuc"))

# This distribution could have been generated from two peaks that were combined. 
# From the histogram, there do not seem to be cells without nuclei or segmentation 
# artifacts where the nucleus is larger than the cell. However, there are so many 
# cells in this dataset and it is possible that just a few cells would not be visible 
# on this histogram. We double check:

# No nucleus
sum(sfe$nucleus_area < 1)
#> [1] 0
# Nucleus larger than cell
sum(sfe$nucleus_area > sfe$cell_area)
#> [1] 0

# So there are no cells without nuclei or nuclei larger than their cells. 
# Here we plot the nuclei proportion in space:
make_plot_grid(fun_voyager = plotSpatialFeature, 
               name_pdf = "prop_nuc_inSpace.pdf",
               main = "nuclei proportion in space",
               list_arg = list(feature = "prop_nuc",
                               colGeometryName = "cellSeg",
                               color = NA))

# Cells in some histological regions have larger proportions occupied by the nuclei. 
# It is interesting to check, controlling for cell type, how cell area, nucleus area, 
# and the proportion of cell occupied by nucleus relate to gene expression. However, 
# a problem in performing such an analysis is that cell segmentation is only available 
# for one z-plane here and these areas also relate to where this z-plane intersects each cell.

# Below we plot a 2D histogram to better show the density of points on this plot:
make_plot_grid(fun_voyager = plotColData, 
               name_pdf = "propNuc_cellArea_scatter.pdf",
               main = "cell area vs nuclei proportion",
               list_arg = list(x = "cell_area",
                               y = "prop_nuc"),
               if_sfe = FALSE,
               name_arg_first = "object")

# Smaller cells tend to have higher proportion occupied by the nucleus. 
# This can be related to cell type, or it could be a limitation in how small
# the nuclei can be in this tissue.

# We also examine the relationship between nucleus area and the proportion of 
# cell occupied by the nucleus:

# plotColData(sfe, x="nucleus_area", y="prop_nuc", bins = 100)
make_plot_grid(fun_voyager = plotColData, 
               name_pdf = "propNuc_nucArea_scatter.pdf",
               main = "nucleus area vs nuclei proportion",
               list_arg = list(x = "nucleus_area",
                               y = "prop_nuc"),
               if_sfe = FALSE,
               name_arg_first = "object")

# The outlier is obvious. There are more cells with both small nuclei and low 
# proportion of area occupied by the nucleus.
# Negative controls
# Since there are only a few hundred genes plus negative control probes, all 
# row names of the SFE object can be printed out to find what the negative control probes are called.

# rownames(sfe)
# According to the Xenium paper (Janesick et al. 2022), there are 3 types of controls:

# probe controls to assess non-specific binding to RNA,
# decoding controls to assess misassigned genes, and
# genomic DNA (gDNA) controls to ensure the signal is from RNA.

# The paper does not explain in detail how those control probes were designed,
# nor explain what the blank probes are. But the blank probes can be used as a 
# negative control.

# add QC-related colData to sfe object for each sample
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
  # colnames(sfe) <- seq_len(ncol(sfe))
  colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area
  # save the updated sfe file with new colData column added
  
  
  is_blank <- str_detect(rownames(sfe), "^BLANK_")
  sum(is_blank)
  # 107 # for sample 6
  
  # This should be number 1, the probe control
  is_neg <- str_detect(rownames(sfe), "^NegControlProbe")
  sum(is_neg)
  #> [1] 20
  
  # This should be number 2, the decoding control
  is_neg2 <- str_detect(rownames(sfe), "^NegControlCodeword")
  sum(is_neg2)
  #> [1] 41
  
  # This must be number 3, gDNA control
  is_anti <- str_detect(rownames(sfe), "^antisense")
  sum(is_anti)
  #> [1] 0
  
  # Also make an indicator of whether a feature is any sort of negative control
  
  is_any_neg <- is_blank | is_neg | is_neg2 | is_anti
  
  
  # The addPerCellQCMetrics() function in the scuttle package can conveniently add transcript counts,
  # proportion of total counts, and number of features detected for any subset of features
  # to the SCE object. Here we do this for the SFE object, as SFE inherits from SCE.
  sfe <- addPerCellQCMetrics(sfe, subsets = list(blank = is_blank,
                                                 negProbe = is_neg,
                                                 negCodeword = is_neg2,
                                                 anti = is_anti,
                                                 any_neg = is_any_neg))
  
  saveRDS(sfe, file = paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
}


# Next we plot the proportion of transcript counts coming from any negative control.
cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")] # get the names of colData columns that are related to the negtive control probes
make_plot_grid(fun_voyager = plotColDataHistogram, 
               name_pdf = "propTrans_negativeControl_histogram.pdf",
               main = "the proportion of transcript counts coming from any negative control",
               list_arg = list(feature = cols_use,
                               bins = 100,
                               ncol=3))
# 7: Removed 33 rows containing non-finite values (`stat_bin()`).
# 8: Removed 275 rows containing non-finite values (`stat_bin()`).
# 9: Removed 22 rows containing non-finite values (`stat_bin()`).
# 10: Removed 66 rows containing non-finite values (`stat_bin()`).
# 11: Removed 11 rows containing non-finite values (`stat_bin()`).
# 12: Removed 22 rows containing non-finite values (`stat_bin()`).

# The histogram is dominated by the bin at zero and there are some extreme outliers
# too few to be seen but evident from the scale of the x axis. We also plot the histogram only 
# for cells with at least 1 count from a negative control. 
# The NA’s come from cells that got segmented but have no transcripts detected.
# plotColDataHistogram(sfe, cols_use, bins = 100, ncol = 3) + 
# scale_x_log10() +
# annotation_logticks(sides = "b")
name_pdf = "propTrans_negativeControl_histogram_log.pdf"
plist_ = make_plot_grid_plist(fun_voyager = plotColDataHistogram, 
                              main = "the proportion of transcript counts coming from any negative control",
                              list_arg = list(feature = cols_use,
                                              bins = 100,
                                              ncol=3))  
for(ksample in 1:length(list_sample)){
  plist_[[ksample]] =  plist_[[ksample]]+
    scale_x_log10() +
    annotation_logticks(sides = "b")
}
pdf(name_pdf, width = 30, height = 20)
print(plot_grid(plotlist = plist_, nrow=2))
dev.off()

# he vast majority of these cells have less than 1% of transcript counts from negative controls, 
# but there are outliers with up to 50%.
# Next we plot the distribution of the number of negative control counts per cell:
cols_use2 <- names(colData(sfe))[str_detect(names(colData(sfe)), "_detected$")]
# plotColDataHistogram(sfe, cols_use2, bins = 20, ncol = 3) +
# Avoid decimal breaks on x axis unless there're too few breaks
# scale_x_continuous(breaks = scales::breaks_extended(Q = c(1,2,5)))
name_pdf = "numTrans_negativeControl_histogram.pdf"
plist_ = make_plot_grid_plist(fun_voyager = plotColDataHistogram, 
                              main = "the number of transcript counts coming from any negative control",
                              list_arg = list(feature = cols_use2,
                                              bins = 20,
                                              ncol=3))  
for(ksample in 1:length(list_sample)){
  plist_[[ksample]] =  plist_[[ksample]]+
    scale_x_continuous(breaks = scales::breaks_extended(Q = c(1,2,5)))   
}
pdf(name_pdf, width = 30, height = 20)
print(plot_grid(plotlist = plist_, nrow=2))
dev.off()
# 
# The counts are low, mostly zero, but there are outliers with up to 10 counts of all types aggregated.
# Then the outlier with 50% of counts from negative controls must have very low total real transcript counts to begin with.
# 
# The scuttle package can detect outliers, but by default it assigns anything above zero as an outlier, 
# since that is over 3 median absolute deviations (MADs) away from the median, which is 0, and the MAD is 0 
# since the vast majority of cells don’t have any negative control count. But it makes sense to allow a small proportion of negative controls. Here we use the distribution just for cells with at least 1 negative control count to find outliers. This distribution has a very long tail and some definite outliers.
# 
# The code below extracts the outliers, based only on cells with at least one negative control count
get_neg_ctrl_outliers <- function(col, sfe) {
  inds <- colData(sfe)$nCounts > 0 & colData(sfe)[[col]] > 0
  df <- colData(sfe)[inds,]
  outlier_inds <- isOutlier(df[[col]], type = "higher")
  outliers <- rownames(df)[outlier_inds]
  col2 <- str_remove(col, "^subsets_")
  col2 <- str_remove(col2, "_percent$")
  new_colname <- paste("is", col2, "outlier", sep = "_")
  colData(sfe)[[new_colname]] <- colnames(sfe) %in% outliers
  sfe
}

# add "is_blank_outlier" to the colData of sfe object for each sample
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
  colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area
  # save the updated sfe file with new colData column added
  
  cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")]
  print(cols_use)
  for (n in cols_use) {
    sfe <- get_neg_ctrl_outliers(n, sfe)
  }
  
  is_blank <- str_detect(rownames(sfe), "^BLANK_")
  sum(is_blank)
  
  sfe <- addPerCellQCMetrics(sfe, subsets = list(is_blank = is_blank ))
  
  saveRDS(sfe, file = paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
}
names(colData(sfe))

# Below we examine where the outliers are located in space:
# plotSpatialFeature(sfe, "is_blank_outlier", colGeometryName = "cellSeg")
make_plot_grid(fun_voyager = plotSpatialFeature, 
               name_pdf = "outlier_negativeControl_InSpace.pdf",
               main = "examine where the outliers are located in space",
               list_arg = list(features = "is_blank_outlier",
                               colGeometryName = "cellSeg",
                               color = NA))

# We find that the outliers are difficult to see:
# plotColData(sfe, y = "is_blank_outlier", x = "cell_area", 
# point_fun = function(...) list()) 
make_plot_grid(fun_voyager = plotColData, 
               name_pdf = "outlierNegativeControl_cell_area_violin.pdf",
               main = "is_blank_outlier vs.cell_area ",
               list_arg = list(y = "is_blank_outlier",
                               x = "cell_area",
                               point_fun = function(...) list()),
               if_sfe = FALSE,
               name_arg_first = "object")

# The analysis reveals that the outliers seem to be smaller. Outliers for negative probe controls 
# and negative codeword controls are also hard to see on the plot, so their plots are skipped here. 
# But the top left region in the tissue tends to have more counts from antisense controls.
# plotSpatialFeature(sfe, "is_anti_outlier", colGeometryName = "cellSeg")
# Yi does not plot this

# Now that we have identified the outliers, we can remove them along with empty cells before proceeding to further analysis:
## generate save sfe_keep (only cells that has passed the QC) for each sample
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_oneSample_",sample_name,".rds"))
  inds_keep <- sfe$nCounts > 0 & sfe$nucleus_area < 400 & !sfe$is_anti_outlier &
    !sfe$is_blank_outlier & !sfe$is_negCodeword_outlier & !sfe$is_negProbe_outlier
  (sfe <- sfe[,inds_keep])
  saveRDS(sfe, file = paste0(dir_data_processed_downstream, "obj_sfe_keep_oneSample_",sample_name,".rds"))
}
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  
  sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_keep_oneSample_",sample_name,".rds"))

  cell_id_kp = colData(sfe)$cell_id
  saveRDS(cell_id_kp, file = paste0(dir_data_processed_downstream, "cell_id_kp_oneSample_",sample_name,".rds"))
  
}
################################################################################################
################################################################################################
################################################################################################
################################################################################################


# Next we check how many negative control features are detected per cell:
# plotColDataHistogram(sfe, cols_use2, bins = 20, ncol = 3) +
# Avoid decimal breaks on x axis unless there're too few breaks
# scale_x_continuous(breaks = scales::breaks_extended(3, Q = c(1,2,5)))
name_pdf = "numberNegativeControl_perCell_histogram.pdf"
plist_ = make_plot_grid_plist(fun_voyager = plotColDataHistogram, 
                              main = "examine how many negative control features are detected per cell",
                              list_arg = list(feature = cols_use2,
                                              bins = 20,
                                              ncol =3),
                              sfe_keep = TRUE)

for(ksample in 1:length(list_sample)){
  plist_[[ksample]] =  plist_[[ksample]]+
    scale_x_continuous(breaks = scales::breaks_extended(Q = c(1,2,5)))   
}
pdf(name_pdf, width = 30, height = 20)
print(plot_grid(plotlist = plist_, nrow=2))
dev.off()


# number of outliers are not as small as it is in the tutorial data
# which means we may need to set stricter cutoff for outlier removal

####################################################################
####################################################################
################################## Genes
####################################################################
####################################################################

# Here we look at the mean and variance of each gene


# add mean, variance and is_neg to the colData for each sample, sfe_keep object
for(ksample in 1:length(list_sample)){
  sample_name = list_sample[ksample]
  sfe = readRDS(paste0(dir_data_processed_downstream, "obj_sfe_keep_oneSample_",sample_name,".rds"))
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  rowData(sfe)$is_neg <- is_any_neg
  saveRDS(sfe, file = paste0(dir_data_processed_downstream, "obj_sfe_keep_oneSample_",sample_name,".rds"))
}
# Real genes generally have higher mean expression across cells than negative controls.

# plotRowData(sfe, x = "means", y = "is_neg") +
# scale_y_log10() +
# annotation_logticks(sides = "b")

name_pdf = "QCgene_mean_vs_negativeControl_violin.pdf"
plist_ = make_plot_grid_plist(fun_voyager = plotRowData, 
                              main = "Real genes generally have higher mean expression across cells than negative controls",
                              list_arg = list(x = "means",
                                              y = "is_neg"),
                              if_sfe = FALSE, 
                              name_arg_first =  "object",
                              sfe_keep = TRUE)

for(ksample in 1:length(list_sample)){
  plist_[[ksample]] =  plist_[[ksample]]+
    scale_y_log10() +
    annotation_logticks(sides = "b")
}
pdf(name_pdf, width = 30, height = 20)
print(plot_grid(plotlist = plist_, nrow=2))
dev.off()


# Here the real genes and negative controls are plotted in different colors
name_pdf = "QCgene_mean_vs_var_scatter_colBy_IfNegativeControl.pdf"
plist_ = make_plot_grid_plist(fun_voyager = plotRowData, 
                              main = "mean vs. variance of each gene",
                              list_arg = list(x = "means",
                                              y = "vars",
                                              color_by = "is_neg"),
                              if_sfe = FALSE, 
                              name_arg_first =  "object",
                              sfe_keep = TRUE)

for(ksample in 1:length(list_sample)){
  plist_[[ksample]] =  plist_[[ksample]]+
    geom_abline(slope = 1, intercept = 0, color = "red") +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() +
    coord_equal() +
    labs(color = "Negative control")
  
}
pdf(name_pdf, width = 30, height = 20)
print(plot_grid(plotlist = plist_, nrow=2))
dev.off()
