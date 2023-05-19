
setwd('/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/')
suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("lobstr")
    library("sessioninfo")
})

## Define some info for the samples
load(here::here("code", "REDCap", "REDCap_HYP.rda"))

# sample_info <- data.frame(dateImg = as.Date(REDCap_dACC$date)) 
# sample_info$experimenterImg <- as.factor(REDCap_dACC$experimenter_img)
sample_info <- data.frame(slide = as.factor(REDCap_HYP$slide))
sample_info$array <- as.factor(REDCap_HYP$array)
sample_info$brnum <- as.factor(sapply(strsplit(REDCap_HYP$sample, "-"), `[`, 1))
sample_info$species <- as.factor(REDCap_HYP$species)
sample_info$replicate <- as.factor(REDCap_HYP$serial)
sample_info$sample_id <- paste(sample_info$slide, sample_info$array, sep = "_")
sample_info$sample_path = file.path(here::here("processed-data", "01_spaceranger", "spaceranger_230511"), sample_info$sample_id,"outs")

##discard barnyard samples
#stopifnot(all(file.exists(sample_info$sample_path)))

#sample_info_human = sample_info[which(sample_info$species == "human"),]
#sample_info_mouse = sample_info[which(sample_info$species == "mouse"),]
# Define the donor info using information from
# donor_info <- read.csv(file.path(here::here("raw-data", "sample_info_Visium", "demographicInfo_Geo.csv")), header = TRUE, stringsAsFactors = FALSE)
# 
# ## check if all donor info included
# setdiff(sample_info_human$brnum,donor_info$brnum)
# 
# ## Combine sample info with the donor info
# sample_info_human <- merge(sample_info_human, donor_info)

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    reference_gtf = file.path("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/","genes", "genes.gtf")
)
Sys.time()
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))

# Sys.time()
# spe <- read10xVisiumWrapper(
#     sample_info_mouse$sample_path,
#     sample_info_mouse$sample_id,
#     type = "sparse",
#     data = "raw",
#     images = c("lowres", "hires", "detected", "aligned"),
#     load = TRUE,
#     reference_gtf = file.path("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-mm10-2020-A/","genes", "genes.gtf")
# )
# Sys.time()
# save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw_human.Rdata"))

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))

## Size in Gb
lobstr::obj_size(spe)
# 5.141452 B
dim(spe)
# 36601 159744

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
