###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
###########
aname <- "normcounts"
lambda=0
###########
dir_data_processed_banksy_old = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/ywang/processed_data/1stRun/"
dir_data_processed_banksy_new = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/ywang/processed_data/2ndRUn/"

###########
dir_out = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/ywang/processed_data/joint/"
setwd(dir_out)
###########
sample_names_old = readRDS(paste0(dir_data_processed_banksy_old, "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))

#############
#############
se = readRDS("spe_joint_2runs_banksy_lambda0_nucleus.rds")

#############
sample_names= c(sample_names_old, sample_names_new)

#############
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )

######################################################################### 
###################################################################### 
## do the plots for cell type clustering results - on 2-d space, and on UMAP
######################################################################### 
library(scater)
plt <- function(se_, 
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  cnames <- colnames(colData(se_))
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  plot_bank1 <- plotColData(se_[,grp1],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank2 <- plotColData(se_[,grp2],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  
  plot_bank_byClu1 <- plotColData(se_[,grp1],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu2 <- plotColData(se_[,grp2],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
    scale_shape_manual(values=seq(0,15))
  
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
  rdnames <- reducedDimNames(se_)
  
  umap_1 <- plotReducedDim(se_[,grp1],
                           dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                           colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 2)+   
    scale_shape_manual(values=seq(0,15))
  
  umap_2<- plotReducedDim(se_[,grp2],
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 2)+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)
  
  print(plot_bank_byClu1)
  print(plot_bank_byClu2)
  
  print(umap_1)
  print(umap_2)
}

pdf("multisample_banksy_2runs_13samples_lambda_0_nuclei.pdf", width = 60, height = 100/2)
# png(paste0("multisample_banksy_2runs_13samples_lambda_0_nuclei.png"), 
    # width = 10, height = 50, units = "in", res = 400)
plt(se)
dev.off()


plt_umap_perClu <- function(se_, 
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  cnames <- colnames(colData(se_))
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  rdnames <- reducedDimNames(se_)
  
  umap_byclu_1 <- plotReducedDim(se_[,grp1],
                                 dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                                 colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))
  
  umap_byclu_2<- plotReducedDim(se_[,grp2],
                                dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                                colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))

  print(umap_byclu_1)
  print(umap_byclu_2)
  
}

# png(paste0(dir_tmp, "/multisample_banksy_2runs_13samples_lambda_0_nuclei_UMAP_perCluster.png"), 
    # width = 10, height = 50, units = "in", res = 400)
pdf("multisample_banksy_2runs_13samples_lambda_0_nuclei_UMAP_perCluster.pdf", width = 60, height = 3*100/2)
plt_umap_perClu(se)
dev.off()

###################################################################### 
## do the plots for gene expression level - on 2-d space
######################################################################### 
#### Specify the genes you want to plot
genes_to_plot = c("ESR1", "GAD1")
########### 

sex_perSampleID = c( c("M", "F", "F", "M", "M", "F"),
                     c("M", "M", "F", "F", "M", "M", "F"))
names(sex_perSampleID) = c(c("br6197", "br5993a", "br5993b", "br1735a", "br1735b", "br6588"),
                           c("br1225a", "br1225b", "br8741c", "br8741d", "br5459a", "br5459b", "br8667c") )
                                                

########### 

# pdf("spatialPlot_nucleus_2runs_together.pdf",width=30,height=30)
for(gene_tmp in genes_to_plot ){
  png(paste0( "spatialPlot_nucleus_2runs_together_gene_", gene, ".png"), 
      width = 10, height = 10, units = "in", res = 400)
  df <- as.data.frame(
    cbind(spatialCoords(se), 
          expr = logcounts(se)[gene_tmp, ]))
  df$sample_id=colData(se)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  
  p <- ggplot(df, aes(x = x_location, y = y_location, color = expr)) + 
    geom_point(size = .0005) + 
    coord_fixed() + 
    scale_color_gradient(low = "gray90", high = "blue", 
                         trans = "sqrt", breaks = range(df$expr), 
                         name = "expression leve") + 
    ggtitle(gene_tmp) + 
    theme_bw() + 
    theme(plot.title = element_text(face = "italic"), 
          # panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())+
    facet_wrap(~sex + sample_id, ncol = 3)+ theme_dark()
  
  print(p)
  dev.off()
}



###################################################################### 
## do the plots for gene expression level - on UMAP
######################################################################### 
library(cowplot)

####
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se[,colData(se)$sizeFactor>0])
se_seurat  <- as.Seurat(se) 
umap_embeddings <- se_seurat@reductions$UMAP_M0_lam0[[1:ncol(se_seurat)]]

#### Specify the genes you want to plot
genes_to_plot = c("ESR1", "GAD1")
  
for(gene in genes_to_plot){
    print(gene)
    png(paste0( "UMAP_nucleus_2runs_together_gene_", gene, ".png"), 
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
    
    
    # Combine the plots into a single figure using cowplot
    combined_plot <- plot_grid(plotlist = plot_list, ncol = 1)
    
    # Display the combined plot
    print(combined_plot)
    
    
    dev.off()
}
  
  




