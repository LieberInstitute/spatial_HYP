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
geneInfo_ = read.csv(paste0(dir_data_processed,"data/Xenium_HYP_Candidate_Genes.csv"))
head(geneInfo_)

###########
########### find genes that are annotated
###########

geneInfo_[,2]

genes_VMH = geneInfo_[grep("VMH",geneInfo_[,2]),]
genes_VMH_marker = genes_VMH[grep("DEG: VMH vs not",genes_VMH[,2]),]
genes_VMH_marker[,2]

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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
se = spe_joint

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
cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]

########################################################################
####### VMH DEG
library(ggpubr)
spe = se
# spe = se[,cellClu=="4"]
# cellClu_used = "4"
if_vmh = ifelse(cellClu=="4", "VMH", "not VMH")

sele = 1:ncol(se)
# sele=sample(1:ncol(se), 10000)
spe = se[,sele]
# spe = 
pdf("AnnotatedVmhMarker_spatialPlot_boxPlot.pdf",width=30,height=30)
for(gene_tmp in intersect(genes_VMH_marker[,1] , rownames(spe))){
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  
  df$sample_id=colData(spe)$sample_id
  df$vmh = if_vmh[sele]
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
  
  df$vmh = factor(df$vmh, levels = c("VMH", "not VMH"))
  # df[df$sex!="Male",]
  p_box <- ggplot(df, aes(x = sex_sampleID, y = expr, fill =  vmh)) + 
    ggtitle(gene_tmp) + 
    theme_bw() + 
    geom_violin(width=1,position=position_dodge(1)) +
    geom_boxplot(width=0.1, position=position_dodge(1)) +

    theme(plot.title = element_text(face = "italic"),
          panel.grid = element_blank()
          # axis.title = element_blank(), 
          # axis.text = element_blank(), 
          # axis.ticks = element_blank()
          )#+
    # scale_fill_manual(values= rep(c( "#56B4E9", "#56B4E9","#E69F00","#E69F00", "#56B4E9","#E69F00"), each = 2))
    # facet_wrap(~sex, ncol = 6)#+ theme_dark()
  
  # print(p_box)
  print(ggarrange(p, p_box, heights = c(2, 0.7),
                  ncol = 1, nrow = 2))
}
dev.off()


