dir_Data_arcPaper = "/dcs04/hansen/data/ywang/xenium/data_paperARC/"

dir_out= "//dcs04/hansen/data/ywang/xenium//integrate_paperARCscRNA/"
setwd(dir_out)

########## load seurat object for scRNA data in ARC paper
obj_arcPaper = readRDS(paste0(dir_Data_arcPaper, "obj_arcPaper_sub.rds"))
# rownames(obj_arcPaper)
dim(obj_arcPaper)

########## load seurat object of xenium data
obj_xenium = readRDS(paste0( "obj_xenium_nuclei_sub.rds"))

dim(obj_xenium)
# [1]    366 10000


########## 
library(Seurat)
# rownames(obj_arcPaper)
# rownames(obj_xenium)

bm280k.list = list(obj_xenium, obj_arcPaper)
features <- SelectIntegrationFeatures(object.list = bm280k.list)


bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = bm280k.list, 
                                  reduction = "rpca",
                                  dims = 1:5
)

bm280k.integrated <- IntegrateData(anchorset = anchors, dims = 1:5)
saveRDS(bm280k.integrated, file = "integrated_subARCpaper_subXemium.rds")

#######

