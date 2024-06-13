# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaName="0"
res=2
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)
library(dplyr)

aname <- "normcounts"


########### load multi-sample Banksy result
se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", res, ".rds"))
# saveRDS(spe_joint, file = paste0("spe_joint_2runs_banksy_lambda0_nucleus.rds"))

# se = spe_joint

################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se[,colData(se)$sizeFactor>0])
# summary(colData(se)$sizeFactor)
se_seurat  <- as.Seurat(se) 
# head(se_seurat@meta.data)
# colData(se)

######################################################################### 
### load unsupervised marker gene result
######################################################################### 

se_seurat_markers = readRDS(paste0(
  "markerGenes_multiSampleBanksy_13samples_lambda",
  lambdaName, "_res",res,".rds"))

se_seurat_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


cellClu_markers = list()
cellClu_uni = unique(se_seurat_markers$cluster)
for(celltype in cellClu_uni){
  
  cellClu_markers[[celltype]] = top10[top10$cluster == celltype,"gene"][[1]]
}
######################################################################### 
### plot the results - UMAP 
######################################################################### 
# Load required libraries
library(Seurat)
library(ggplot2)


umap_embeddings <- se_seurat@reductions$UMAP_M0_lam0[[1:ncol(se_seurat)]]

for(cellType in cellClu_uni[1:20]){
  dir_tmp = paste0("unsupervisedMarkerGenes_umap_clu", cellType)
  dir.create(dir_tmp)
  print(cellType)
  
  # Specify the genes you want to plot
  genes_to_plot = cellClu_markers[[cellType]]
  
  print(cellType)
  
  print(genes_to_plot)
  for(gene in genes_to_plot){
    png(paste0(dir_tmp, "/", gene, ".png"), 
        width = 10, height = 10, units = "in", res = 400)

    # Loop through each gene and create a plot
    plot_list <- lapply(gene, function(gene) {
      # Create a data frame with UMAP coordinates, cell identities, and gene expression
      plot_data <- FetchData(se_seurat, vars = c("sample_id", gene))
      # plot_data$ident <- Idents(se_seurat)
      plot_data$UMAP_1 = umap_embeddings[,1]
      plot_data$UMAP_2 = umap_embeddings[,2]
      
      # Create the plot using ggplot2
      plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = get(gene))) +
        geom_point(size = 0.5) +
        scale_color_gradient(low = "gray90", high = "red") +
        facet_wrap(~sample_id, ncol = 3) +  # Split plots by cell identity
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        labs(title = gene, x = "UMAP 1", y = "UMAP 2", color = "Expression")
      
      return(plot)
    })
    # print(p)
    
    
    # Combine the plots into a single figure using cowplot
    library(cowplot)
    combined_plot <- plot_grid(plotlist = plot_list, ncol = 1)
    
    # Display the combined plot
    print(combined_plot)
  # }
  
    dev.off()
    
  }
}

######################################################################### 
######################################################################### 
######################################################################### 
for(cellType in cellClu_uni[1:41]){
  print(cellType)
  dir_tmp = paste0("unsupervisedClu_umap")
  dir.create(dir_tmp)

    png(paste0(dir_tmp, "/clu", cellType, ".png"), 
        width = 10, height = 10, units = "in", res = 400)
    
    # Loop through each gene and create a plot
    # plot_list <- lapply(gene, function(gene) {
      # Create a data frame with UMAP coordinates, cell identities, and gene expression
      plot_data <- FetchData(se_seurat, vars = c("sample_id", "clust_M0_lam0_k50_res2"))
      # plot_data$ident <- Idents(se_seurat)
      plot_data$UMAP_1 = umap_embeddings[,1]
      plot_data$UMAP_2 = umap_embeddings[,2]
      plot_data$cluster = plot_data[,"clust_M0_lam0_k50_res2"]
      plot_data$within_clu = as.character( plot_data$cluster == cellType)
      
      # Create the plot using ggplot2
      plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = within_clu)) +
        geom_point(size = 0.5) +
        # scale_color_gradient(low = "gray90", high = "red") +
        facet_wrap(~sample_id, ncol = 3) +  # Split plots by cell identity
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        labs(title = paste0("cluster ", cellType), x = "UMAP 1", y = "UMAP 2")

    print(plot)

    dev.off()
    
  #
}
