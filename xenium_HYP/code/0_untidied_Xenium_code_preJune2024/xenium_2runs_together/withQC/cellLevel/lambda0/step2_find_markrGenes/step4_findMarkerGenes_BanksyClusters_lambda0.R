######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R


lambdaName=0
lambda=0
res <- 2

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed_banksy =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)
###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
# library(Banksy)
aname <- "normcounts"

###########
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")

se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", res, ".rds"))

######################################################################### 
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se[,colData(se)$sizeFactor>0])
se <- scater::logNormCounts(se)
se_seurat  <- as.Seurat(se) 

######################################################################### 
### identify markergenes for each cl
######################################################################### 

cnames <- colnames(colData(se))

cnames <- cnames[grep("^clust", cnames)]
print(cnames)
  
  Idents(se_seurat) = colData(se)[,cnames]
  
  #####
  se_seurat_markers = FindAllMarkers(se_seurat)
  
  saveRDS(se_seurat_markers, file = paste0(
    "markerGenes_multiSampleBanksy_13samples_lambda",
    lambdaName, "_res",res,".rds"))
  ######################################################################### 
  se_seurat_markers = readRDS(paste0(
    "markerGenes_multiSampleBanksy_13samples_lambda",
    lambdaName, "_res",res,".rds"))
  ######################################################################### 
  ### plot the results - heatmap for marker genes
  ######################################################################### 
  ### genes to check for cl 33
  library(dplyr)
  genes_interest = c("gad1", "gad2", "slc32a1", "slc17a6",  "slc17a7",
                     "KISS1","TAC3","GHRH","GSX1","SST","PCP4","AGRP",
                     "NPY","POMC","CBLN1","TBX3")
  genes_interest = toupper(genes_interest)
  genes_interest = intersect(genes_interest, rownames(se))
  # > genes_interest
  # [1] "GAD1"    "GAD2"    "SLC17A6" "SLC17A7" "TAC3"    "GHRH"    "GSX1"   
  # [8] "SST"     "PCP4"    "AGRP"    "POMC"  
  
  ######## format se into seurat format
  obj = CreateSeuratObject(assays(se)$counts, project = paste0("marker genes"))
  # obj@data = 
  obj = NormalizeData(obj)
  # obj@data 
  obj = ScaleData(obj)
  Idents(obj) = colData_[,paste0("clust_M0_lam", lambda, "_k50_res",res)]
  
  ######## 
  clus_unique = unique( Idents(se_seurat) )
  nclu = length(unique( Idents(se_seurat) ))
  # Idents(se_seurat) = factor(Idents(se_seurat) , levels= as.character(c(33,1:32,34:nclu)))
  
  ######## get top 10 marker genes for each cell cluster
  se_seurat_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  ######## rd sample cells
  set.seed(111) #plot a subset of cells to reduce running time
  rd_cells = sample(1:ncol(obj), 20000)
  # features_plot = c(genes_interest,(top10$gene)[!(top10$gene)%in%genes_interest])
  
  ######## create pdf
  # patchwork_1.1.3
  # library(patchwork)
  library(pheatmap)
  library(grid)
  
  clus_unique = as.character(clus_unique[order(as.integer(clus_unique))])
  
  for(k in 1:length(clus_unique)){
    clu = clus_unique[k]
    print(clu)

    # print(DoHeatmap(obj[features_plot, rd_cells] , features =features_plot) + NoLegend())
    features_plot_tmp = top10[top10$cluster== as.character(clu) ,"gene"][[1]]
    Idents(se_seurat) = factor(Idents(se_seurat) , 
                               levels=  c(clu, clus_unique[-k]) )
    # 
    if(clu == "33"){features_plot_tmp = c(genes_interest, features_plot_tmp ) }
    
    mat_tmp = data.matrix(obj$RNA@scale.data)[as.character(features_plot_tmp), rd_cells]
    
    col_annotation = data.frame(cluster = Idents(se_seurat) [rd_cells])
    rownames(col_annotation) = colnames(obj)[rd_cells]

    
    p = pheatmap(mat_tmp[,order(col_annotation$cluster)], show_colnames=F, main=paste0("cluster ", clu),
                 annotation_col = col_annotation,cluster_cols=F,cluster_rows=F)
    
    pdf(paste0("markerGenes_top10_multisampleBanksy_lambda_",
               lambdaName,
               "_res",
               res,"_clu",clu,".pdf"), width = 20, height = 20)
    print(p)
    dev.off()
    
  }

  ############################ scatter plot
  mat_tmp = data.matrix(obj$RNA@scale.data)[c(toupper("slc17a6"), "GAD1"), Idents(se_seurat)=="33"]
  pdf("scatter_slc17a6_gad1_cl33.pdf")
  slc17a6 = mat_tmp[toupper("slc17a6"),]
  GAD1 = mat_tmp[toupper("GAD1"),]
  plot(GAD1,slc17a6)
  dev.off()

  
  ############################ top 50
  ############################
  se_seurat_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 50) %>%
    ungroup() -> top50
  
  write.csv(top50, file=paste0(
    "top50_markerGenes_multiSampleBanksy_13samples_lambda",
    lambdaName,
    "_res", res, ".csv"))

  