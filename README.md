# spatial_HYP


##Internal: `/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/`


## /Code Directories:

###01_Spaceranger: Initial data handling through Spaceranger (Ryan M)

###02_build_spe: Initial assembly of spatialExperiment object (Ryan M)

###03_filtering: QC and filtering on technical variables (mitochondrial reads, number genes detected, etc) (Bernie M)


###REDCAP: ?

###Vistoseg: Align section image data with spaceranger (Ryan M) 



## /analysis/code directories:
00-Run Harmony.rmd: Correcting rownames of colData and colnames of spe to be unique for downstream/future BioC compatability; run Harmony by brnum (since technical replicates are present from some donors).

###01-feature_selection: HVGs, nnSVG, others(? future)

###02-clustering: walktrap based clustering on sets of top x percentile HVGs, nnSVGs, and combinations thereof, across SNNgraph k of 10, 15...30

###H01-feature_selection: feature selection as above, on Harmony-corrected data.

###H02-clustering: clustering as above but with walktrap k=5, 10, 15 only, on Harmony-corrected data. 