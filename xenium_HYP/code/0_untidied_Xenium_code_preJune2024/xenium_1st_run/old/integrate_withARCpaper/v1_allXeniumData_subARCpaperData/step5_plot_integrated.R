dir_Data_arcPaper = "/dcs04/hansen/data/ywang/xenium/data_paperARC/"

dir_out= "//dcs04/hansen/data/ywang/xenium//integrate_paperARCscRNA/"
setwd(dir_out)

####### load integrated data

bm280k.integrated  = readRDS("integrated_subARCpaper.rds")

####### calculate PCA for integrated data

library(Seurat)
bm280k.integrated = ScaleData(bm280k.integrated)
bm280k.integrated = FindVariableFeatures(bm280k.integrated)
bm280k.integrated = RunPCA(bm280k.integrated, npcs=5)

bm280k.integrated = readRDS("integrated_subARCpaper_withPCA.rds")
bm280k.integrated = RunUMAP(bm280k.integrated, dims=1:5)
saveRDS(bm280k.integrated, file = "integrated_subARCpaper_withPCA_UMAP.rds")


####### plot integrated data
bm280k.integrated = readRDS("integrated_subARCpaper_withPCA_UMAP.rds")

dim(bm280k.integrated)
# [1]    349 376888

DimPlot(bm280k.integrated, 
        split.by="orig.ident",
        # plot.title="integrated_hm_scRNA_adr_subcells_v2.rds",
        reduction="pca",ncol=2)

bm280k.integrated@meta.data$if_clu15  = as.character(Idents(bm280k.integrated)=="15") 
DimPlot(bm280k.integrated[,], 
        split.by =c("if_clu15"),
        reduction="pca",ncol=3)


# plot umap

DimPlot(bm280k.integrated, 
        split.by="orig.ident",
        # plot.title="integrated_hm_scRNA_adr_subcells_v2.rds",
        reduction="umap",ncol=2)

DimPlot(bm280k.integrated[,], 
        split.by =c("if_clu15"),
        reduction="umap",ncol=2)
