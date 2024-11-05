suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

spe_vis_in <- here(
    "processed-data", "13-SamuiSPE.RDS"
)

spe_out <- here(
    "processed-data", "14-samui", "spe.h5ad"
)


spe_r_out <- here(
    "processed-data", "14-samui", "spe_r.rds"
)

spe <- readRDS(spe_vis_in)

# load(spe_vis_in)

###############################################################################
#  Functions
###############################################################################

write_anndata <- function(spe, out_path) {
    invisible(
        basiliskRun(
            fun = function(spe, filename) {
                library("zellkonverter")
                library("reticulate")
                
                # Convert SCE to AnnData:
                adata <- SCE2AnnData(spe)
                
                #  Write AnnData object to disk
                adata$write(filename = filename)
                
                return()
            },
            env = zellkonverterAnnDataEnv(),
            spe = spe,
            filename = out_path
        )
    )
}

###############################################################################
#   Main
###############################################################################

#-------------------------------------------------------------------------------
#   Convert snRNA-seq and spatial R objects to AnnData python objects
#-------------------------------------------------------------------------------

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R
reducedDims(spe)$spatial <- spatialCoords(spe)

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(spe) <- rowData(spe)$gene_id

#   Save a copy of the filtered + slightly modified sce as an R object, and
#   convert all objects to Anndatas
# saveRDS(spe, spe_r_out)

print("Converting objects to AnnDatas...")
write_anndata(spe, spe_out)

