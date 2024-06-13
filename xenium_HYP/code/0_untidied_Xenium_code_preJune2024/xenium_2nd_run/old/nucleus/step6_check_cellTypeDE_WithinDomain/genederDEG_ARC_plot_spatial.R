# ml conda_R/4.3
# R
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
library(ggpubr)
library(ggplot2)
library(cowplot)
########### 
aname <- "normcounts"

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("M", "F", "F", "M", "M", "F")
names(sex_perSampleID) = sampleID_unique

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
setwd(dir_out)
###################################################################### 
### load stuff
###################################################################### 
########### load cell type clustering using nucle onlu
cellType_clu_nuc = readRDS("cellType_clu_nuc.rds")
# saveRDS(cellType_clu_nuc, file = "cellType_clu_nuc.rds")
###########

###########
geneInfo_ = read.csv(paste0(dir_data_processed,"data/Xenium_HYP_Candidate_Genes.csv"))

########### load cell-type-domain level DE significancy result
# write.csv(mat_bin_overall, file="mat_mat_bin_overall_perCellType_ARC_overall.csv",quote=F)
mat_bin_overall = read.csv("mat_mat_bin_overall_perCellType_ARC_overall.csv")
head(mat_bin_overall)
###########
########### find genes that are annotated
###########


genes_ARC = geneInfo_[grep("ARC",geneInfo_[,2]),]
genes_ARC_sexDEG = genes_ARC[grep("Sex-differential",genes_ARC[,2]),]


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
                       cellIDs_se_domain[cellClu_domain == "3"])
# length(cellIDs_kp)

########### 
cellClu_domain_kp = colData(se_domain)[cellIDs_kp,paste0("clust_M0_lam", 1, "_k50_res0.6")]
cellClu_annotation_kp = colData(se)[cellIDs_kp,paste0("clust_M0_lam", 0, "_k50_res0.6")]



########################################################################
################################# format data
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]



########################################################################
########################################################################
####### ARC DEG

spe = se[,cellIDs_kp]
cellType_clu_unique= sort(unique(cellClu_annotation_kp))


# genes_ARC_sexDEG_kp = intersect(genes_ARC_sexDEG[-1,1],
                                # rownames(spe))
# length(genes_ARC_sexDEG_kp)
# 14
cellTypeAnnotations = c("astrocytes", "likely more than one cell types","myelinating oligodendrocytes",
                        "microglia","","oligodendrocyte precursor",
                        "immune related","epithelial cells","astrocytes",
                        "excitatory neurons","", "POMC neuronal cells",
                        "non-myelinating matureish(?) oligodendrocytes","epithelial cells","",
                        "","T cells","AVP neurons",
                        "monocyte")
names(cellTypeAnnotations) = as.character(1:19)
# > length(genes_ARC_sexDEG_kp)
# [1] 14


######## 
# SHANK3
######## plot, folder by cell gene

library(ggpubr)

for(gene_tmp in genes_ARC_sexDEG[,1] ){
  print(gene_tmp)
  
  plist = list()
  
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  
  dir_geneTmp = paste0("ARC_spatial_plot_",gene_tmp)
  dir.create(dir_geneTmp)
  # dir.create(paste0("spatial_plot_",gene_tmp))
  for(cellType_clu in cellType_clu_unique){
    print(cellType_clu)
    pdf(paste0(dir_geneTmp,"/ARC_DEG_spatialperCellClu",cellType_clu,".pdf"),width=30,height=30)
    

    df_cellType_clu_tmp = df[cellClu_annotation_kp == cellType_clu,]
    
    
    # plist[[as.integer(cellType_clu)]]   
    p <- ggplot(df_cellType_clu_tmp, aes(x = x_location, y = y_location, color = expr)) +
    geom_point(size = .0005) +
    coord_fixed() +
    scale_color_gradient(low = "gray90", high = "blue",
                         trans = "sqrt", breaks = range(df$expr),
                         name = "expression level") +
      ggtitle(paste0(gene_tmp,", cluster ", cellType_clu,", ", cellTypeAnnotations[as.character(cellType_clu)])) + 
      theme_bw() +
    theme(plot.title = element_text(face = "italic"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    facet_wrap(~sex + sample_id, ncol = 3)+ theme_dark()
    
    
    p_box <- ggplot(df_cellType_clu_tmp, aes(x = sex, y = expr, fill = sample_id)) + 
      # geom_boxplot(position=position_dodge(1))+
      # geom_point(size = .0005) +
      # coord_fixed() +
      # scale_color_gradient(low = "gray90", high = "blue", 
      #                      trans = "sqrt", breaks = range(df$expr), 
      #                      name = "expression leve") + 
      ggtitle(gene_tmp) + 
      theme_bw() + 
      geom_violin(width=1,position=position_dodge(1)) +
      geom_boxplot(width=0.1, position=position_dodge(1)) +
      # geom_jitter()+
      # shape=16, position=position_jitter(0.2)
      theme(plot.title = element_text(face = "italic"),
            panel.grid = element_blank()
            # axis.title = element_blank(), 
            # axis.text = element_blank(), 
            # axis.ticks = element_blank()
      )+
      scale_fill_manual(values=c( "#56B4E9", "#56B4E9","#E69F00","#E69F00", "#56B4E9","#E69F00"))
    
    
    print(ggarrange(p, p_box, heights = c(2, 0.7),
                    ncol = 1, nrow = 2))
    
    dev.off()
    
  }
  # print(plot_grid(plotlist=plist, ncol=4))
}

######## plot, folder by cell cluster
library(ggpubr)

for(gene_tmp in genes_ARC_sexDEG[,1] ){
  print(gene_tmp)
  
  plist = list()
  
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  

  # dir.create(paste0("spatial_plot_",gene_tmp))
  for(cellType_clu in cellType_clu_unique){
    
    dir_clusterTmp = paste0("ARC_spatial_plot_cluster",cellType_clu)
    dir.create(dir_clusterTmp)
    
    print(cellType_clu)
    pdf(paste0(dir_clusterTmp,"/ARC_DEG_spatial_",gene_tmp,".pdf"),width=30,height=30)
    
    
    df_cellType_clu_tmp = df[cellClu_annotation_kp == cellType_clu,]
    
    
    # plist[[as.integer(cellType_clu)]]   
    p <- ggplot(df_cellType_clu_tmp, aes(x = x_location, y = y_location, color = expr)) +
      geom_point(size = .0005) +
      coord_fixed() +
      scale_color_gradient(low = "gray90", high = "blue",
                           trans = "sqrt", breaks = range(df$expr),
                           name = "expression level") +
      ggtitle(paste0(gene_tmp,", cluster ", cellType_clu,", ", cellTypeAnnotations[as.character(cellType_clu)])) + 
      theme_bw() +
      theme(plot.title = element_text(face = "italic"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())+
      facet_wrap(~sex + sample_id, ncol = 3)+ theme_dark()
    
    
    p_box <- ggplot(df_cellType_clu_tmp, aes(x = sex, y = expr, fill = sample_id)) + 
      # geom_boxplot(position=position_dodge(1))+
      # geom_point(size = .0005) +
      # coord_fixed() +
      # scale_color_gradient(low = "gray90", high = "blue", 
      #                      trans = "sqrt", breaks = range(df$expr), 
      #                      name = "expression leve") + 
      ggtitle(gene_tmp) + 
      theme_bw() + 
      geom_violin(width=1,position=position_dodge(1)) +
      geom_boxplot(width=0.1, position=position_dodge(1)) +
      # geom_jitter()+
      # shape=16, position=position_jitter(0.2)
      theme(plot.title = element_text(face = "italic"),
            panel.grid = element_blank()
            # axis.title = element_blank(), 
            # axis.text = element_blank(), 
            # axis.ticks = element_blank()
      )+
      scale_fill_manual(values=c( "#56B4E9", "#56B4E9","#E69F00","#E69F00", "#56B4E9","#E69F00"))
    
    
    print(ggarrange(p, p_box, heights = c(2, 0.7),
                    ncol = 1, nrow = 2))
    
    dev.off()
    
  }
  # print(plot_grid(plotlist=plist, ncol=4))
}
