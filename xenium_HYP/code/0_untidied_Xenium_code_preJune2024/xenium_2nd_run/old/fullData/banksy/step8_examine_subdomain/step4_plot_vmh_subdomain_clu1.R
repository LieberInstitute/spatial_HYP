# ml conda_R/4.3
# R

###########
lambda <- 1
lambdaName="1"

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
setwd(dir_out)


###########
library(ggpubr)
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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_vmhOnly.rds"))
se = spe_joint

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("M", "F", "F", "M", "M", "F")
names(sex_perSampleID) = sampleID_unique




########################################################################
########################################################################
# saveRDS(cluster1.markers, 
        # file= "cluster1_markers_vmh.rds")
cluster1.markers = readRDS("cluster1_markers_vmh.rds")
se_seurat_markers = cluster1.markers
se_seurat_markers %>%
  # group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  # slice_head(n = 10) %>%
  ungroup() -> top10

########################################################################
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)
colData(se)[,which(colnames(colData(se))=="clust_M0_lam1_k50_res0.6")[1]] = NULL # rm the layer annotations using full data

cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]
spe = se
########################################################################

####### VMH subdoman 1 marker gene
# dir.create("VMH_subdomain1")
# pdf("VMH_subdomain1_markerGenes_spatialPlot_boxplot.pdf",width=30,height=30,units="in)
clu_tmp="1"
for(gene_tmp in rownames(top10) ){
  
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  df$subdomain = as.character(cellClu)
  df$subdomain[df$subdomain != clu_tmp] = paste0("outside subdomain ", clu_tmp) 
  df$subdomain[df$subdomain == clu_tmp] = paste0("subdomain ", clu_tmp)
  
  df$subdomain = factor(df$subdomain, levels = c(paste0("subdomain ", clu_tmp), 
                                                 paste0("outside subdomain ", clu_tmp) ))
  
  png(paste0("VMH_subdomain1/",gene_tmp,".png"),res=200,
        width=30,height=30,units="in")
  # png(paste0(dir_out_clu_tmp, "/", gene_tmp, ".png"), res=200,
      # width=30, height=30, units="in")
  
  
  range_ = range(df$expr)
  xlim_ = range(df$x_location)
  ylim_ = range(df$y_location)
  
  p2 <- ggplot(df[df$subdomain != paste0("subdomain ", clu_tmp), ], aes(x = x_location, y = y_location, color = expr)) +
    geom_point(size = .0005) +
    coord_fixed() +
    scale_color_gradient(low = "white", high = "blue",
                         trans = "sqrt",
                         breaks = range_, 
                         limits = range_,
                         name = "expression level") +
    ggtitle(paste0(gene_tmp,", outside subdomain ", clu_tmp) )+
    # theme_bw() +
    xlim(xlim_)+ylim(ylim_)+
    theme(plot.title = element_text(size =55), 
          panel.background = element_rect(fill = "gray24",
                                          colour = "gray24"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    facet_wrap(~sex + sample_id, ncol = 3)
  
  
  
  p <- ggplot(df[df$subdomain == paste0("subdomain ", clu_tmp), ], aes(x = x_location, y = y_location, color = expr)) + 
    geom_point(size = .0005) + 
    coord_fixed() + 
    xlim(xlim_)+ylim(ylim_)+
    scale_color_gradient(low = "white", high = "blue", 
                         trans = "sqrt", 
                         breaks = range_, 
                         limits = range_,
                         name = "expression level") + 
    ggtitle(paste0(gene_tmp,", inside subdomain ", clu_tmp) )+
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
  
  p_box <- ggplot(df, aes(x = subdomain, y = expr, fill = subdomain)) + 
    
    ggtitle(gene_tmp) + 
    theme_bw() + 
    geom_violin(width=1,position=position_dodge(1)) +
    geom_boxplot(width=0.1, position=position_dodge(1)) +
    theme(plot.title = element_text(face = "italic"),
          panel.grid = element_blank()  )+
    facet_wrap(~sex + sample_id, ncol = 3)#+
  
  p_box2 <- ggplot(df, aes(x = subdomain, y = expr, fill = subdomain)) + 
    
    ggtitle(gene_tmp) + 
    theme_bw() + 
    geom_violin(width=1,position=position_dodge(1)) +
    geom_boxplot(width=0.1, position=position_dodge(1)) +
    theme(plot.title = element_text(face = "italic"),
          panel.grid = element_blank()  )+
    facet_wrap(~sex + sample_id, ncol = 6)#+
  
  print(ggarrange(p, p2,p_box2, heights = c(2, 2, .7),
                  ncol = 1, nrow = 3))
  dev.off()
  
  
}

# for(gene_tmp in rownames(top10) ){
#   
#   png(paste0("VMH_subdomain1/",gene_tmp,".png"),res=200,
#       width=30,height=30,units="in")
#   
#   df <- as.data.frame(
#     cbind(spatialCoords(spe), 
#           expr = logcounts(spe)[gene_tmp, ]))
#   df$sample_id=colData(spe)$sample_id
#   
#   df$sex=sex_perSampleID[as.character(df$sample_id)]
#   df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
#   df$subdomain = as.character(cellClu)
#   df$subdomain[df$subdomain !="1"] = "outside subdomain 1" 
#   df$subdomain[df$subdomain =="1"] = "subdomain 1" 
#   
#   df$subdomain = factor(df$subdomain, levels = c("subdomain 1", "outside subdomain 1"))
#   p <- ggplot(df, aes(x = x_location, y = y_location, color = expr)) + 
#     geom_point(size = .0005) + 
#     coord_fixed() + 
#     scale_color_gradient(low = "gray90", high = "blue", 
#                          trans = "sqrt", breaks = range(df$expr), 
#                          name = "expression leve") + 
#     ggtitle(gene_tmp) + 
#     theme_bw() + 
#     theme(plot.title = element_text(face = "italic"), 
#           # panel.grid = element_blank(), 
#           axis.title = element_blank(), 
#           axis.text = element_blank(), 
#           axis.ticks = element_blank())+
#     facet_wrap(~sex + sample_id, ncol = 3)+ theme_dark()
#   
#   p_box <- ggplot(df, aes(x = subdomain, y = expr, fill = subdomain)) + 
# 
#     ggtitle(gene_tmp) + 
#     theme_bw() + 
#     geom_violin(width=1,position=position_dodge(1)) +
#     geom_boxplot(width=0.1, position=position_dodge(1)) +
#     # geom_jitter()+
#     # shape=16, position=position_jitter(0.2)
#     theme(plot.title = element_text(face = "italic"),
#           panel.grid = element_blank()  )+
#     facet_wrap(~sex + sample_id, ncol = 3)#+
#     # scale_fill_manual(values=c( "#56B4E9", "#56B4E9","#E69F00","#E69F00", "#56B4E9","#E69F00"))
#   # +
#   # facet_wrap(~sex, ncol = 6)#+ theme_dark()
#   print(ggarrange(p, p_box, heights = c(2, 0.7),
#                   ncol = 1, nrow = 2))
#   # print(p+p_box)
#   dev.off()
#   
#   
# }



