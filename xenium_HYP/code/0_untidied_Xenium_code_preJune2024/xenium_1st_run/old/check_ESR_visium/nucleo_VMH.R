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
# geneInfo_ = read.csv(paste0(dir_data_processed,"data/Xenium_HYP_Candidate_Genes.csv"))
# head(geneInfo_)

###########
########### find genes that are annotated
###########

# geneInfo_[,2]
# 
# genes_VMH = geneInfo_[grep("VMH",geneInfo_[,2]),]
# genes_VMH_sexDEG = genes_VMH[grep("sex",genes_VMH[,2]),]
# genes_VMH_sexDEG[,2]

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

########### load multi-sample Banksy result, cell-level
lambda <- 1
lambdaName="1"
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
se = spe_joint
cellClu_domain = colData(spe_joint)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
se_domain = spe_joint
cellIDs_se_domain = colData(se_domain)$cell_id
names(cellClu_domain) = cellIDs_se_domain
########### 

############ load cell type clustering result, nucleus only
lambda <- 0
lambdaName="0"
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_nucleus.rds"))
se = spe_joint
cellID_cellClustAnnotation=colData(se)$cell_id



###################################################################### 
########### keep cells of interest
###################################################################### 
########### keep cells that are in VMH domain and are kept in both process - full data and nucleo
cellIDs_kp = intersect(cellID_cellClustAnnotation,
                       cellIDs_se_domain[cellClu_domain=="4"])

########### 
cellClu_domain_kp = colData(se_domain)[cellIDs_kp,paste0("clust_M0_lam", 1, "_k50_res0.6")]
cellClu_annotation_kp = colData(se)[cellIDs_kp,paste0("clust_M0_lam", 0, "_k50_res0.6")]

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female")
names(sex_perSampleID) = sampleID_unique

########### 
# sampleID_unique <- c("br6197","br5993a","br5993b",
#               "br1735a","br1735b","br6588")
# sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female")
# names(sex_perSampleID) = sampleID_unique

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
# spe = se
spe = se[,cellIDs_kp]
# cellClu_used = "4"


# "SNHG14" %in% rownames(spe) - some genes may have beem filtered out in the QC step
# FALSE

pdf("ESR1_VMH_nucleus.pdf",width=30,height=30)
for(gene_tmp in "ESR1"){
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  
  df$sample_id=colData(spe)$sample_id
  
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
  
  # print(p)
  
  
  
  p_box <- ggplot(df, aes(x = sex, y = expr, fill = sample_id)) + 
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
  # +
    # facet_wrap(~sex, ncol = 6)#+ theme_dark()
  
  print(ggarrange(p, p_box, heights = c(2, 0.7),
                  ncol = 1, nrow = 2))  
}
dev.off()
# sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female")
# col = c("#E69F00", "#56B4E9")
# names(col)=c("Male","Female")




