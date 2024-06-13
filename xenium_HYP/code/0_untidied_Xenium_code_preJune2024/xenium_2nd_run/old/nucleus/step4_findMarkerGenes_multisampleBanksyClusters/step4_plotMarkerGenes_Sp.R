# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaName="0"

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
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
aname <- "normcounts"



###################################################################### 
################################### load multi-sample Banksy result
###################################################################### 
############ load domain segmentation result
lambda <- 1
lambdaName="1"
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
cellClu_domain = colData(spe_joint)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
se_domain = spe_joint
cellIDs_se_domain = colData(se_domain)$cell_id
names(cellClu_domain) = cellIDs_se_domain

############ load cell type clustering result, nucleus only
lambda <- 0
lambdaName="0"
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_nucleus.rds"))
se = spe_joint
cellID_cellClustAnnotation=colData(se)$cell_id


###################################################################### 
########### keep cells of interest
###################################################################### 
########### keep cells that are in ARC domain and are kept in both process - full data and nucleo
cellIDs_kp = intersect(cellID_cellClustAnnotation,
                       cellIDs_se_domain[cellClu_domain=="3"])

########### 
cellClu_domain_kp = colData(se_domain)[cellIDs_kp,paste0("clust_M0_lam", 1, "_k50_res0.6")]
cellClu_annotation_kp = colData(se)[cellIDs_kp,paste0("clust_M0_lam", 0, "_k50_res0.6")]

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female")
names(sex_perSampleID) = sampleID_unique


########################################################################
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)
# cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

# se_seurat  <- as.Seurat(se)
# se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
# sex = sex_perSampleID[colData(se)$sample_id]

########################################################### 
####### ARC DEG
######################################################################### 
library(ggpubr)

spe = se[,cellIDs_kp]


######################################################################### 
################################# load marker genes of each cell cluster
######################################################################### 


library(dplyr)
se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda",lambdaName,"_nucleus.rds"))

se_seurat_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


######## plot, folder by gene
# RNF152
library(ggpubr)
clu_unique = as.character(unique(cellClu_annotation_kp))

for(kClu in 1:length(clu_unique)){
  clu_tmp = clu_unique[kClu]
  print(clu_tmp)
  
  clusterTmp_markers = top10[top10$cluster==clu_tmp,"gene"]$gene
  
  dir_out_clu_tmp = paste0("cellCluster_markers_Space_plotARConly_cluster", clu_tmp)
  dir.create(dir_out_clu_tmp)
  
  for(gene_tmp in clusterTmp_markers){
    
    df <- as.data.frame(
      cbind(spatialCoords(spe), 
            expr = logcounts(spe)[gene_tmp, ]))
    df$sample_id=colData(spe)$sample_id
    
    df$sex=sex_perSampleID[as.character(df$sample_id)]
    df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
    
    df$cellClu = as.character(cellClu_annotation_kp)
    
    df$with_cellClu = ""
    df$with_cellClu[df$cellClu != clu_tmp] = paste0("outside cell cluster ", clu_tmp) 
    df$with_cellClu[df$cellClu == clu_tmp] = paste0("within cell cluster ", clu_tmp)
    
    df$with_cellClu = factor(df$with_cellClu, levels = c(paste0("within cell cluster ", clu_tmp), 
                                                   paste0("outside cell cluster ", clu_tmp) ))
    
    png(paste0(dir_out_clu_tmp, "/", gene_tmp, ".png"), res=200,
        width=30, height=30, units="in")
    
    
    range_ = range(df$expr)
    xlim_ = range(df$x_location)
    ylim_ = range(df$y_location)
    
    p2 <- ggplot(df[df$cellClu != paste0(clu_tmp), ], 
                 aes(x = x_location, y = y_location, color = expr)) +
      geom_point(size = .0005) +
      coord_fixed() +
      scale_color_gradient(low = "white", high = "blue",
                           trans = "sqrt",
                           breaks = range_, 
                           limits = range_,
                           name = "expression level") +
      ggtitle(paste0(gene_tmp,", outside cell cluster ", clu_tmp) )+
      # theme_bw() +
      xlim(xlim_)+ylim(ylim_)+
      theme(plot.title = element_text(size =55), 
            panel.background = element_rect(fill = "gray24",
                                            colour = "gray24"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())+
      facet_wrap(~sex + sample_id, ncol = 3)
    
    
    
    p <- ggplot(df[df$cellClu == paste0(clu_tmp), ], 
                aes(x = x_location, y = y_location, color = expr)) + 
      geom_point(size = .0005) + 
      coord_fixed() + 
      xlim(xlim_)+ylim(ylim_)+
      scale_color_gradient(low = "white", high = "blue", 
                           trans = "sqrt", 
                           breaks = range_, 
                           limits = range_,
                           name = "expression level") + 
      ggtitle(paste0(gene_tmp,", inside cell cluster ", clu_tmp) )+
      # theme_bw() + 
      theme(plot.title = element_text(size =55), 
            panel.background = element_rect(fill = "gray24",
                                            colour = "gray24",
                                            size = 0.5),
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())+
      facet_wrap(~sex + sample_id, ncol = 3)
    # + theme_dark()
    
    p_box <- ggplot(df, aes(x = with_cellClu, y = expr, fill = with_cellClu)) + 
      
      ggtitle(gene_tmp) + 
      theme_bw() + 
      geom_violin(width=1,position=position_dodge(1)) +
      geom_boxplot(width=0.1, position=position_dodge(1)) +
      theme(plot.title = element_text(face = "italic"),
            panel.grid = element_blank()  )+
      facet_wrap(~sex + sample_id, ncol = 3)#+
    
    p_box2 <- ggplot(df, aes(x = with_cellClu, y = expr, fill = with_cellClu)) + 
      
      ggtitle(gene_tmp) + 
      theme_bw() + 
      geom_violin(width=1,position=position_dodge(1)) +
      geom_boxplot(width=0.1, position=position_dodge(1)) +
      theme(plot.title = element_text(face = "italic"),
            panel.grid = element_blank()  )+
      facet_wrap(~sex + sample_id, ncol = 6)#+
    
    print(ggarrange(p, p2, p_box2, heights = c(2, 2, .7),
                    ncol = 1, nrow = 3))
    dev.off()
    
    
  }
}
