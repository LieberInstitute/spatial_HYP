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
cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]
spe = se
########################################################################
# genes_interested =c("gad1", "gad2", "slc17a6", "slc17a7", "slc32a1") 
genes_interested = toupper(genes_interested)

subdomain_unique =unique( as.character(cellClu))
####### VMH subdoman 1 marker genesubdomain_unique =unique( as.character(cellClu))

# dir.create("VMH_GenesOfInterested")
# pdf("VMH_subdomain1_markerGenes_spatialPlot_boxplot.pdf",width=30,height=30,units="in)

genes_interested_kp  = intersect(genes_interested, rownames(spe))
# slc32a1 is not in the xenium data
for(gene_tmp in genes_interested_kp ){

  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  df$subdomain = as.character(cellClu)
  # df$subdomain[df$subdomain !="1"] = "outside subdomain 1" 
  # df$subdomain[df$subdomain =="1"] = "subdomain 1" 
  
  # df$subdomain = factor(df$subdomain, levels = c("subdomain 1", "outside subdomain 1"))
  range_value_exp = range(df$expr)
  
  pdf(paste0("VMH_GenesOfInterested/",gene_tmp,".pdf"),
      # ,res=200,
      width=30,height=30)
  # ,units="in")
  
  for(kSubdomain in 1:length(subdomain_unique)){
    
    subdomain_tmp  = sort(subdomain_unique)[kSubdomain]
    p <- ggplot(df[df$subdomain==subdomain_tmp,], aes(x = x_location, y = y_location, color = expr)) + 
      geom_point(size = .0005) + 
      coord_fixed() + 
      scale_color_gradient(low = "gray90", high = "blue", 
                           trans = "sqrt", breaks = range(df$expr), 
                           name = "expression level",
                           limits=range_value_exp) + 
      ggtitle(paste0(gene_tmp,", subdomain ", subdomain_tmp)) + 
      theme_bw() + 
      theme(plot.title = element_text(face = "italic"), 
            # panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())+
      facet_wrap(~sex + sample_id, ncol = 3)+ theme_dark()
    
    
    print(p)
    
    
  }
  dev.off()
  
  
}



