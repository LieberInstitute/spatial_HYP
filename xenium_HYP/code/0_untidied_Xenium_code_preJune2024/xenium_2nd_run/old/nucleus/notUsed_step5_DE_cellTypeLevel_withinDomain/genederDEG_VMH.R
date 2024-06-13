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
genes_VMH_sexDEG = genes_VMH[grep("sex",genes_VMH[,2]),]
genes_VMH_sexDEG[,2]

# [1] "Sex-differential gene of ARC cluster 8 and clusters 8+13 combined"                                                               
# [2] "Sex-differential gene of ARC cluster 8"                                                                                          
# [3] "Sex-differential gene of ARC cluster 8, possible marker of ARC clusters 8+13"                                                    
# [4] "Sex-differential gene of ARC cluster 8, possible marker of ARC clusters 8+13 (slash of cluster 8 as relativistic marker over 13)"
# [5] "Sex-differential gene of ARC cluster 8 and clusters 8+13 combined"                                                               
# [6] "Sex-differential gene of ARC cluster 8; possible ARC cl8 specific marker?"                                                       
# [7] "Sex-differential gene of ARC cluster 8, mouse ARC cell type marker (POMC neurons), ARC spot marker (this data)"                  
# [8] "Sex-differential gene of ARC cluster 8 and clusters 8+13 combined"                                                               
# [9] "Sex-differential gene of ARC cluster 8; gene with high spot level expression in VMH"                                             
# [10] "Sex-differential gene of ARC cluster 8 and 8+13; sex-differential gene of VMH"                                                   
# [11] "Sex-differential gene of ARC cluster 8; possible spatial pattern"                                                                
# [12] "Sex-differential gene of ARC cluster 8; possible marker (by expression level) for ARC"                                           
# [13] "Sex-differential gene of ARC cluster 8+13 combined; possible spatial pattern"                                                    
# [14] "Sex-differential gene of ARC cluster 8 and clusters 8+13 combined"                                                               
# [15] "Sex-differential gene of ARC cluster 8"         

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
####### ARC DEG
library(ggpubr)

spe = se[,cellClu=="4"]
cellClu_used = "4"


# "SNHG14" %in% rownames(spe) - some genes may have beem filtered out in the QC step
# FALSE

pdf("AnnotatedSexDEG_VMH_spatialPlot_boxPlot.pdf",width=30,height=30)
for(gene_tmp in intersect(genes_VMH_sexDEG[-1,1] , rownames(spe))){
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




