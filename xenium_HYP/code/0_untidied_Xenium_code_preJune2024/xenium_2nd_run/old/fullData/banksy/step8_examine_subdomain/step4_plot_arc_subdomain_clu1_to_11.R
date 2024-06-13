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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_arcOnly.rds"))
se = spe_joint

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("M", "F", "F", "M", "M", "F")
names(sex_perSampleID) = sampleID_unique


########################################################################
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)
colData(se)[,which(colnames(colData(se))=="clust_M0_lam1_k50_res0.6")[1]] = NULL # rm the layer annotations using full data
# Idents(obj) = colData(se)[,"clust_M0_lam1_k50_res0.6"]

cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]
spe = se
########################################################################

####### VMH subdoman 1 marker gene
# dir.create("VMH_subdomain1")
clu_unique = as.character(sort(as.integer(unique(cellClu))))

mat_markerGene = array("", dim=c(10, length(clu_unique)))
colnames(mat_markerGene) = paste0("cluster ", clu_unique)
# mat_markerGeCne
for(kClu in 1:length(clu_unique)){
  clu_tmp = clu_unique[kClu]
  print(clu_tmp)
  
  cluster_tmp.markers = readRDS(paste0("cluster",clu_tmp,"_markers_arc.rds"))
  
  cluster_tmp.markers %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  if(nrow(top10)>0){
    mat_markerGene[1:nrow(top10),paste0("cluster ",clu_tmp)] = rownames(top10)
    
  }
}


write.csv(mat_markerGene, file="markerGene_each_ARC_subdomain.csv", quote=F)
# pdf("VMH_subdomain1_markerGenes_spatialPlot_boxplot.pdf",width=30,height=30,units="in)
for(kClu in 1:length(clu_unique)){
  clu_tmp = clu_unique[kClu]
  print(clu_tmp)
  
  cluster_tmp.markers = readRDS(paste0("cluster",clu_tmp,"_markers_arc.rds"))
  
  cluster_tmp.markers %>%
    # group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  dir_out_clu_tmp = paste0("ARC_subdomain", clu_tmp)
  # dir.create(dir_out_clu_tmp)
  
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
    
    png(paste0(dir_out_clu_tmp, "/", gene_tmp, ".png"), res=200,
        width=30, height=30, units="in")
    

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
}










