# ml conda_R/4.3
# R
###########
res=2
###########
lambda <- 0
lambdaName="0"

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


# Load required libraries
library(Seurat)
library(ggplot2)
cellClu_markers = list()
cellClu_uni = unique(se_seurat_markers$cluster)
for(celltype in cellClu_uni){

    cellClu_markers[[celltype]] = top10[top10$cluster == celltype,"gene"][[1]]
}

######################################################################### 
### plot the results - UMAP 
######################################################################### 

library(cowplot)
# 18 cell types in total
for(cellType in cellClu_uni[36:41]){
  print(cellType)
  genes_to_plot = cellClu_markers[[cellType]]
  
  nrow_ = length(genes_to_plot)/2
  
  dir_tmp = paste0("markerGenes_unsupervised_spacePlot_clu", cellType)
  dir.create(dir_tmp)
  
  print("cellType")
  
  print(cellType)

  for(gene in genes_to_plot){
    png(paste0(dir_tmp, "/", gene, ".png"), 
        width = 10, height = 10, units = "in", res = 200)
    # x = "x_location", y = "y_location",
    
    # Loop through each gene and create a plot
    plot_list <- lapply(gene, function(gene) {
      # Create a data frame with UMAP coordinates, cell identities, and gene expression
      plot_data <- FetchData(se_seurat, vars = c("sample_id", gene))
      # plot_data$ident <- Idents(se_seurat)
      plot_data$x_position = colData(se)$x_position
      plot_data$y_position = colData(se)$y_position
      
      # Create the plot using ggplot2
      plot <- ggplot(plot_data, aes(x = x_position, y = y_position, color = get(gene))) +
        geom_point(size = 0.0003, alpha=0.1) +
        scale_color_gradient(low = "gray90", high = "red") +
        facet_wrap(~sample_id, ncol = 3) +  # Split plots by cell identity
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        labs(title = gene, x = "x_position", y = "y_position", color = "Expression")
      
      return(plot)
    })
    
    
    # Combine the plots into a single figure using cowplot
    library(cowplot)
    combined_plot <- plot_grid(plotlist = plot_list, ncol = 1)
    
    # Display the combined plot
    print(combined_plot)
    # }
    
    dev.off()
    
  }
  

    
}



