# ml conda_R/4.3
# R

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"

###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
dir_data_processed =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
###########
setwd(dir_out)
###########
geneInfo_ = read.csv(paste0(dir_data_processed,"Xenium_HYP_Candidate_Genes.csv"))
head(geneInfo_)

###########
########### find genes that are annotated
###########

genes_ARC = geneInfo_[grep("ARC",geneInfo_[,2]),]
genes_ARC_sexDEG = genes_ARC[grep("Sex-differential",genes_ARC[,2]),]

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
se  = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res2.rds"))

########### 
########### 
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names= c(sample_names_old, sample_names_new)
sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female",
                    c("Male","Male","Female","Female","Male","Male","Female"))
names(sex_perSampleID) = sample_names

########################################################################
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)
cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res2")]

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]


########################################################################
knn_pred_arc_c_ = readRDS("knn_ARC_twoRounds.rds")

########################################################################
####### plot using boxplot
library(ggpubr)
spe = se[,knn_pred_arc_c_=="ARC"]

# > dim(genes_ARC_sexDEG)
# [1] 15  7

pdf("AnnotatedSexDEG_ARC_spatialPlot_boxplot_cellLevel.pdf",width=30,height=30)
for(gene_tmp in genes_ARC_sexDEG[,1] ){
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  
  ### remove one sample where it seems has little ARC contained
  df = df[df$sample_id!="br8667c",]
  
  ### assign levels to the factor denoting sample ID based on sex - male and female (for better organized in the plotting)
  SampleIDs_kp = unique(df$sample_id)
  nsample_kp = length(SampleIDs_kp)
  df$sample_id = factor(df$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
  df$sex = factor(df$sex, level = c("Male", "Female"))
  
  ### do the plot
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
  
  
  p_box <- ggplot(df, aes(x = sex, y = expr, fill = sample_id)) + 
    ggtitle(gene_tmp) + 
    theme_bw() + 
    geom_violin(width=1,position=position_dodge(1)) +
    geom_boxplot(width=0.1, position=position_dodge(1)) +
    theme(plot.title = element_text(face = "italic"),
          panel.grid = element_blank()
          # axis.title = element_blank(), 
          # axis.text = element_blank(), 
          # axis.ticks = element_blank()
    )+
    scale_fill_manual(values=c( "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9","#56B4E9","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00" ))
  print(ggarrange(p, p_box, heights = c(2, 0.7),
                  ncol = 1, nrow = 2))

}
dev.off()





