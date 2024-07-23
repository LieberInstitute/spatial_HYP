suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("data.table")
    library("sessioninfo")
})

## make sure we're in the visium working directory
stopifnot(length(grep(here(),pattern="spatial_HYP/spatial_HYP",value=T))==1)


## Define some info for the samples
## not all data is in redcap, so use a manual table instead.

demos <- fread("raw-data/demos_by_sampleid.txt")
demos[,species:="Human"]
demos[,replicate:=1]
## indicate samples which are duplicates
demos[sample_id %in% paste0("V12D05-350_",c("C1","D1")),replicate:=2]

## we only have one replicate for Br5459--we aren't using the first as it was re-run with better quality results on a second slide.
demos[brnum=="Br5459",replicate:=2]

## extract slide and array from sample id
demos[,slide:=gsub(sample_id,pattern="^(.*)_.*$",replacement="\\1")]
demos[,array:=gsub(sample_id,pattern="^.*_(.*)$",replacement="\\1")]

## first 8 samples are in a different directory than the last couple samples
demos[!(slide %in% c("V13M13-362","V13Y24-346")),sample_path:=here("processed-data/01_spaceranger/spaceranger_230511",sample_id,"outs")]

## more recent samples:
demos[slide %in% c("V13M13-362","V13Y24-346"),sample_path:=here("processed-data/01_spaceranger",sample_id,"outs")]

demos <- data.frame(demos)

## Build basic SPE
Sys.time()
# [1] "2024-07-22 21:41:23 EDT"

spe <- read10xVisiumWrapper(
        demos$sample_path,
    demos$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    reference_gtf = file.path("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/","genes", "genes.gtf")
)
# 2024-07-22 21:41:23.184923 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2024-07-22 21:45:46.129528 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2024-07-22 21:46:02.630822 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2024-07-22 21:46:06.58463 rtracklayer::import: reading the reference GTF file
# 2024-07-22 21:47:32.12815 adding gene information to the SPE object
# 2024-07-22 21:47:32.477257 adding information used by spatialLIBD

Sys.time()
# [1] "2024-07-22 21:47:42 EDT"

## make  colnames unique (spe$key)
colnames(spe) <- spe$key

## append study info to colData, get rid of 10x stuff that was pulled into colData and reducedDims that we don't care about

tmpcd <- as.data.table(colData(spe),keep.rownames=T)
keepnames <- grep(names(tmpcd),pattern="^X10x",value=T,invert=T)
tmpcd <- tmpcd[,..keepnames]
tmpcd <- merge.data.table(tmpcd,as.data.table(demos),by="sample_id",all.x=T)
setnames(tmpcd,c("Best.RIN.PFC","BMI..calculated."),c("best_RIN_PFC","bmi"))
tmpcd <- DataFrame(tmpcd,row.names=tmpcd$rn)[colnames(spe),]
tmpcd$rn <- NULL
colData(spe) <- tmpcd

### 10x dimred garbage
reducedDims(spe) <- reducedDims(spe)[!grepl("10x",names(reducedDims(spe)))]

## Read in cell counts and segmentation results
segmentations_list <-
  lapply(demos$sample_path, function(p) {
    file <- paste0(p,
        "/spatial",
        "/tissue_spot_counts.csv"
      )
    if (!file.exists(file)) {
      return(NULL)
    }
    x <- read.csv(file)
    sampleid <- gsub(p,pattern="^/dcs04/.*/(V.*-.*_.*)/outs$",replacement="\\1")
    x$key <- paste0(x$barcode, "_", sampleid)
    return(x)
  })

## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
  Reduce(function(...) {
    merge(..., all = TRUE)
  }, segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
  segmentations[segmentation_match, -which(
    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
  )]
colData(spe) <- cbind(colData(spe), segmentation_info)

## Redefine tissue overlap from loupe for vis. later with spatialLIBD if desired
spe$overlaps_tissue <-
  factor(ifelse(spe$in_tissue, "in", "out"))

dim(spe)
## raw dimensions: 36601, 49842

## save the raw SPE before removing blank features, nontissue spots, etc
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw_n10.Rdata"))

## now minimally filter to a starting SPE (no undetected genes, tissue spots only)

## drop the spots outside the tissue
spe <- spe[, spe$in_tissue]
dim(spe)
# [1] 36601 45649

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
  message("removing spots without counts for spe")
  spe <- spe[, -which(colSums(counts(spe)) == 0)]
  dim(spe)
}

# removing spots without counts for spe
# [1] 36601 45641

## Remove genes with no data (ie in tissue spots)
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 6238
spe <- spe[-no_expr, ]

dim(spe)
### output dims: 30363, 45641

saveRDS(spe, file = here::here("processed-data", "02_build_spe", "spe_n10.RDS"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2024-07-22 22:15:26 EDT"
# user   system  elapsed
# 949.505   13.538 2537.980
# ─ Session info ──────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 Patched (2024-02-08 r85876)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-07-22
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────
# package                * version     date (UTC) lib source
# abind                    1.4-5       2016-07-21 [2] CRAN (R 4.3.2)
# AnnotationDbi            1.64.1      2023-11-03 [2] Bioconductor
# AnnotationHub            3.10.0      2023-10-24 [2] Bioconductor
# attempt                  0.3.1       2020-05-03 [2] CRAN (R 4.3.2)
# beachmat                 2.18.0      2023-10-24 [2] Bioconductor
# beeswarm                 0.4.0       2021-06-01 [2] CRAN (R 4.3.2)
# benchmarkme              1.0.8       2022-06-12 [2] CRAN (R 4.3.2)
# benchmarkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.2)
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel             1.36.0      2023-10-24 [2] Bioconductor
# BiocSingular             1.18.0      2023-10-24 [2] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
# Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
# cachem                   1.1.0       2024-05-16 [1] CRAN (R 4.3.2)
# cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# codetools                0.2-19      2023-02-01 [3] CRAN (R 4.3.2)
# colorspace               2.1-0       2023-01-23 [2] CRAN (R 4.3.2)
# config                   0.3.2       2023-08-30 [2] CRAN (R 4.3.2)
# cowplot                  1.1.3       2024-01-22 [2] CRAN (R 4.3.2)
# crayon                   1.5.3       2024-06-20 [1] CRAN (R 4.3.2)
# curl                     5.2.1       2024-03-01 [1] CRAN (R 4.3.2)
# data.table             * 1.15.4      2024-03-30 [1] CRAN (R 4.3.2)
# DBI                      1.2.3       2024-06-02 [1] CRAN (R 4.3.2)
# dbplyr                   2.5.0       2024-03-19 [1] CRAN (R 4.3.2)
# DelayedArray             0.28.0      2023-10-24 [2] Bioconductor
# DelayedMatrixStats       1.24.0      2023-10-24 [2] Bioconductor
# digest                   0.6.34      2024-01-11 [2] CRAN (R 4.3.2)
# doParallel               1.0.17      2022-02-07 [2] CRAN (R 4.3.2)
# dotCall64                1.1-1       2023-11-28 [2] CRAN (R 4.3.2)
# dplyr                    1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
# dqrng                    0.4.1       2024-05-28 [1] CRAN (R 4.3.2)
# DropletUtils             1.22.0      2023-10-24 [1] Bioconductor
# DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
# edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
# fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
# fastmap                  1.2.0       2024-05-15 [1] CRAN (R 4.3.2)
# fields                   15.2        2023-08-17 [2] CRAN (R 4.3.2)
# filelock                 1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
# foreach                  1.5.2       2022-02-02 [2] CRAN (R 4.3.2)
# generics                 0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb           * 1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11      2024-02-09 [2] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges          * 1.54.1      2023-10-29 [2] Bioconductor
# ggbeeswarm               0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2                  3.5.1       2024-04-23 [1] CRAN (R 4.3.2)
# ggrepel                  0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
# glue                     1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
# golem                    0.4.1       2023-06-05 [2] CRAN (R 4.3.2)
# gridExtra                2.3         2017-09-09 [2] CRAN (R 4.3.2)
# gtable                   0.3.5       2024-04-22 [1] CRAN (R 4.3.2)
# HDF5Array                1.30.0      2023-10-24 [2] Bioconductor
# here                   * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# htmltools                0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets              1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                   1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
# httr                     1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
# interactiveDisplayBase   1.40.0      2023-10-24 [2] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [2] Bioconductor
# irlba                    2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
# iterators                1.0.14      2022-02-05 [2] CRAN (R 4.3.2)
# jquerylib                0.1.4       2021-04-26 [2] CRAN (R 4.3.2)
# jsonlite                 1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
# KEGGREST                 1.42.0      2023-10-24 [2] Bioconductor
# later                    1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
# lattice                  0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval                 0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
# lifecycle                1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
# limma                    3.58.1      2023-10-31 [2] Bioconductor
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# magick                   2.8.3       2024-02-18 [1] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
# Matrix                   1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
# matrixStats            * 1.3.0       2024-04-11 [1] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# mime                     0.12        2021-09-28 [2] CRAN (R 4.3.2)
# munsell                  0.5.1       2024-04-01 [1] CRAN (R 4.3.2)
# paletteer                1.6.0       2024-01-21 [2] CRAN (R 4.3.2)
# pillar                   1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig                2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
# plotly                   4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
# png                      0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
# promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# purrr                    1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R.methodsS3              1.8.2       2022-06-13 [2] CRAN (R 4.3.2)
# R.oo                     1.26.0      2024-01-24 [2] CRAN (R 4.3.2)
# R.utils                  2.12.3      2023-11-18 [2] CRAN (R 4.3.2)
# R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                     1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
# rhdf5                    2.46.1      2023-11-29 [2] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters             1.14.1      2023-11-06 [2] Bioconductor
# Rhdf5lib                 1.24.2      2024-02-07 [2] Bioconductor 3.18 (R 4.3.2)
# rjson                    0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
# rlang                    1.1.4       2024-06-04 [1] CRAN (R 4.3.2)
# rprojroot                2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools                2.18.0      2023-10-24 [2] Bioconductor
# RSQLite                  2.3.7       2024-05-27 [1] CRAN (R 4.3.2)
# rstudioapi               0.15.0      2023-07-07 [2] CRAN (R 4.3.2)
# rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer            * 1.62.0      2023-10-24 [2] Bioconductor
# S4Arrays                 1.2.1       2024-03-04 [1] Bioconductor 3.18 (R 4.3.2)
# S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
# scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scater                   1.30.1      2023-11-16 [2] Bioconductor
# scuttle                  1.12.0      2023-10-24 [2] Bioconductor
# sessioninfo            * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
# shiny                    1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
# shinyWidgets             0.8.1       2024-01-10 [2] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0      2023-10-24 [2] Bioconductor
# spam                     2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray              1.2.4       2024-02-11 [1] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats        1.14.0      2023-10-24 [2] Bioconductor
# SpatialExperiment      * 1.12.0      2023-10-24 [2] Bioconductor
# spatialLIBD            * 1.14.1      2023-11-30 [2] Bioconductor 3.18 (R 4.3.2)
# statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
# tibble                   3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                    1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect               1.2.1       2024-03-11 [1] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
# vctrs                    0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
# vipor                    0.4.7       2023-12-18 [2] CRAN (R 4.3.2)
# viridis                  0.6.5       2024-01-29 [2] CRAN (R 4.3.2)
# viridisLite              0.4.2       2023-05-02 [2] CRAN (R 4.3.2)
# XML                      3.99-0.16.1 2024-01-22 [2] CRAN (R 4.3.2)
# xtable                   1.8-4       2019-04-21 [2] CRAN (R 4.3.2)
# XVector                  0.42.0      2023-10-24 [2] Bioconductor
# yaml                     2.3.8       2023-12-11 [2] CRAN (R 4.3.2)
# zlibbioc                 1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.2)
#
# [1] /users/bmulvey/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
# [4] /users/bmulvey/R/4.3
#
# ─────────────────────────────────────────────────────────────────────────────
