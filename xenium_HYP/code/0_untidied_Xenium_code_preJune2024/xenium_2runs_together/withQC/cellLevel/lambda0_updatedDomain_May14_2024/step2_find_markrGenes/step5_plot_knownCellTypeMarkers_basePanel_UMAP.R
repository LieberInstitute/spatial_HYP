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

table_basePanel = read.csv("/dcs04/hansen/data/ywang/xenium_joint_runs/Xenium_hBrain_v1_metadata.csv")
# Genes      Ensembl_ID Num_Probesets Codewords Annotation
# 1    ABCC9 ENSG00000069431             8         1       VLMC
# 2 ADAMTS12 ENSG00000151388             8         1       VLMC
# 3 ADAMTS16 ENSG00000145536             8         1      L4 IT
# 4  ADAMTS3 ENSG00000156140             8         1      L6 IT
# 5   ADRA1A ENSG00000120907             8         1       Sncg
# 6   ADRA1B ENSG00000170214             8         1       Sncg
dim(table_basePanel)
# > dim(table_basePanel)
# [1] 266   5

table(table_basePanel$Annotation)
cellType_uni = unique(table_basePanel$Annotation)
cellType_uni = cellType_uni[!cellType_uni%in%c("Glioblastoma (Cancer cells)", "L2/3 IT", "L4 IT", "L5 ET", "L6 CT", "L6 IT", "L6 IT Car3",
                                               "L6b", "L5 IT", "L5/6 NP")]

sum(table_basePanel$Annotation %in% cellType_uni)
# 198
table_basePanel = table_basePanel[table_basePanel$Annotation %in% cellType_uni,]
# Genes
# Annotation
# table_basePanel
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
celltype_markers = list()
for(celltype in cellType_uni){
  which_ = which(table_basePanel$Annotation==celltype)
  celltype_markers[[celltype]] = intersect(table_basePanel[which_,"Genes"],
                                               rownames(se))
}



celltype_markers = remove_empty_items(celltype_markers)


## code below: claude 3
umap_embeddings <- se_seurat@reductions$UMAP_M0_lam0[[1:ncol(se_seurat)]]


library(cowplot)
# 18 cell types in total
for(cellType in cellType_uni[13]){
  print(cellType)
  genes_to_plot = celltype_markers[[cellType]]
  
  nrow_ = length(genes_to_plot)/2
  
  dir_tmp = paste0("markerGenes_umap_", cellType)
  dir.create(dir_tmp)
  
  print("cellType")
  
  print(cellType)
  print("genes_to_plot")
  
  print(genes_to_plot)
  for(gene in genes_to_plot){
    print(gene)
    png(paste0(dir_tmp, "/UMAP_nucleus_basePanel_2runs_together_gene_", gene, ".png"), 
        width = 10, height = 10, units = "in", res = 400)
    # Specify the genes you want to plot
    
    
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
    
    
    # Combine the plots into a single figure using cowplot
    combined_plot <- plot_grid(plotlist = plot_list, ncol = 1)
    
    # Display the combined plot
    print(combined_plot)
    
    
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
  return(filtered_list)
}

# Assuming you have a Seurat object named 'seurat_object' with the UMAP coordinates and cell identities

# genes_to_plot <- c("OXT", "AVP", "CRH", "POMC", "AGRP")




