# ml conda_R/4.3
# R

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
se  = readRDS(paste0("/dcs04/hansen/data/ywang/xenium_joint_runs/out/spe_joint_2runs_banksy_lambda0_nucleus.rds"))
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
### plot the results - UMAP 
######################################################################### 
# Load required libraries
library(Seurat)
library(ggplot2)
hypothalamus_markers <- list(
  magnocellular_neurosecretory_cells = c("OXT", "AVP"),
  parvocellular_neurosecretory_cells = c("CRH", "TRH", "GNRH1", 
                                         "SST", "TH", "SLC6A3"),
  NPY_neurons = "NPY",
  AgRP_neurons = "AGRP",
  POMC_neurons = "POMC",
  kisspeptin_neurons = "KISS1",
  orexin_neurons = "HCRT",
  MCH_neurons = "PMCH",
  galanin_neurons = "GAL",
  histaminergic_neurons = "HDC",
  glutamatergic_neurons = c("SLC17A6", "SLC17A7"),
  GABAergic_neurons = c("GAD1", "GAD2", "SLC32A1"),
  cholinergic_neurons = "CHAT",
  serotonergic_neurons = c("TPH2", "SLC6A4"),
  NOS_expressing_neurons = "NOS1",
  # glial_cells = list(
    astrocytes = c("GFAP", "S100B", "AQP4"),
    oligodendrocytes = c("OLIG2", "MBP", "CNP"),
    microglia = c("CX3CR1", "IBA1", "ITGAM"),
    tanycytes = c("RAXD2", "DIO2", "CNTFR")
  # )
)
for(celltype in names(hypothalamus_markers)){
  hypothalamus_markers[[celltype]] = intersect(hypothalamus_markers[[celltype]],
                                               rownames(se))
}



hypothalamus_markers = remove_empty_items(hypothalamus_markers)


## code below: claude 3
# umap_embeddings <- se_seurat@reductions$UMAP_M0_lam0[[1:ncol(se_seurat)]]
# umap_embeddings = unlist(umap_embeddings)
# head(x = FetchData(object = se_seurat, vars.fetch = 'UMAP1'))
# FetchData(object = se_seurat, vars.fetch = 'UMAP')
# cellType = names(hypothalamus_markers)[1]

for(cellType in names(hypothalamus_markers)){
  dir_tmp = paste0("markerGenes_spacePlot_", cellType)
  dir.create(dir_tmp)
  
  # Specify the genes you want to plot
  genes_to_plot = hypothalamus_markers[[cellType]]
  
  print(cellType)
  
  print(genes_to_plot)
  for(gene in genes_to_plot){
    png(paste0(dir_tmp, "/spacePlot_nucleus_basePanel_2runs_together_gene_", gene, ".png"), 
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
        geom_point(size = 0.0003,alpha=0.1) +
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

  
  
#################################################
################################################# 
# below is from claude 
remove_empty_items <- function(lst) {
  # Use lapply() to iterate over the list elements
  filtered_list <- lapply(lst, function(item) {
    # Check if the item is a list itself
    if (is.list(item)) {
      # Recursively call the function on nested lists
      item <- remove_empty_items(item)
    }
    # Return the item if it is not empty or NULL
    if (!is.null(item) && length(item) > 0) {
      return(item)
    }
  })
  
  # Remove NULL elements from the filtered list
  filtered_list <- filtered_list[!sapply(filtered_list, is.null)]
  
  # Return the filtered list
}

# Assuming you have a Seurat object named 'seurat_object' with the UMAP coordinates and cell identities

# genes_to_plot <- c("OXT", "AVP", "CRH", "POMC", "AGRP")




