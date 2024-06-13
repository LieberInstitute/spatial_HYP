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
lambdaName=0

aname <- "normcounts"
########### 
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names= c(sample_names_old, sample_names_new)
sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female",
                    c("Male","Male","Female","Female","Male","Male","Female"))
names(sex_perSampleID) = sample_names

###########
# dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
# dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
# dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
# dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out_nucleus/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")
setwd(dir_out)
dir_data_processed =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"

###################################################################### 
### load stuff
###################################################################### 
########### load cell type clustering using nucle onlu
# cellType_clu_nuc = readRDS("cellType_clu_nuc.rds")
# saveRDS(cellType_clu_nuc, file = "cellType_clu_nuc.rds")
###########
# write.csv(mat_adj_bin_overall_withFC, file="mat_adj_bin_overall_withFC_ARC.csv",quote=F)
mat_adj_bin_overall_withFC = read.csv("mat_adj_bin_overall_withFC_ARC.csv")
###########
# geneInfo_ = read.csv(paste0(dir_data_processed,"/Xenium_HYP_Candidate_Genes.csv"))

########### load cell-type-domain level DE significancy result
# write.csv(mat_bin_overall, file="mat_mat_bin_overall_perCellType_ARC_overall.csv",quote=F)
# mat_bin_overall = read.csv("mat_mat_bin_overall_perCellType_ARC_overall.csv")
# head(mat_bin_overall)
###########
########### find genes that are annotated
###########


# genes_ARC = geneInfo_[grep("ARC",geneInfo_[,2]),]
# genes_ARC_sexDEG = genes_ARC[grep("Sex-differential",genes_ARC[,2]),]
# 

###################################################################### 
################################### load multi-sample Banksy result
###################################################################### 

############ load cell type clustering result, nucleus only
lambda <- 0
lambdaName="0"
se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", 2, ".rds"))

cellClu_annotation = colData(se)[,"clust_M0_lam0_k50_res2"]
length(cellClu_annotation)
# 879396
########### get ARC domain
spe_joint=se
if_ARC = array("not ARC", dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("/dcs04/hansen/data/ywang/xenium_joint_runs/out/knn_predARC_",sample,"_k500_2ndRun.rds"))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
  if_ARC[which_sample][knn_pred>0.3] = "ARC"
}

length(if_ARC)
# 879396


###################################################################### 
########### keep cells of interest
###################################################################### 
########### keep cells that are in ARC domain and are kept in both process - full data and nucleo
# cellIDs_kp = intersect(cellID_cellClustAnnotation,
                       # cellIDs_se_domain[cellClu_domain == "3"])
# # length(cellIDs_kp)
# 
# ########### 
# # cellClu_domain_kp = colData(se_domain)[cellIDs_kp,paste0("clust_M0_lam", 1, "_k50_res0.6")]
# cellClu_annotation_kp = colData(se)[cellIDs_kp,paste0("clust_M0_lam", 0, "_k50_res0.6")]
# 


########################################################################
################################# format data
######################################################################### 
######################################################################### 
se <- scater::computeLibraryFactors(se)
##
which_cells_kp =  which(colData(se)$sizeFactor>0 & if_ARC =="ARC")

##
cellClu_domain_kp = if_ARC[which_cells_kp]
cellClu_annotation_kp =cellClu_annotation[which_cells_kp]
##
se <- scater::logNormCounts(se[,which_cells_kp])

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]


########################################################################
########################################################################
####### ARC DEG

# spe = se[,cellIDs_kp]
cellType_clu_unique= sort(unique(cellClu_annotation_kp))
nclu = length(cellType_clu_unique)

# genes_ARC_sexDEG_kp = intersect(genes_ARC_sexDEG[-1,1],
                                # rownames(spe))
# length(genes_ARC_sexDEG_kp)
# 14
# cellTypeAnnotations = c("astrocytes", "likely more than one cell types","myelinating oligodendrocytes",
#                         "microglia","","oligodendrocyte precursor",
#                         "immune related","epithelial cells","astrocytes",
#                         "excitatory neurons","", "POMC neuronal cells",
#                         "non-myelinating matureish(?) oligodendrocytes","epithelial cells","",
#                         "","T cells","AVP neurons",
#                         "monocyte")
cellTypeAnnotations = 1:nclu
names(cellTypeAnnotations) = as.character(1:nclu)
# > length(genes_ARC_sexDEG_kp)
# [1] 14


######## 
# SHANK3
######## plot, folder by cell gene
# > genes_ARC_sexDEG[,1]
# [1] "ADARB1"   "ALKAL2"   "SHISAL2B" "CRH"      "DCC"      "IL18R1"  
# [7] "PROKR1"   "RDH10"    "RXRG"     "SHANK3"   "SSTR2"    "VGF"     
# [13] "HRH1"     "LYPD6"    "C2orf80" 
spe = se


spe_cl30 = se[,cellClu_annotation_kp==30]
sex_perCell_cl30 = sex_perSampleID[colData(spe_cl30)$sample_id]
colData(spe_cl30)
which_f= sex_perCell_cl30=="Female"
which_m= sex_perCell_cl30=="Male"

# spe$data
logcounts_spe_cl30 = assays(spe_cl30)$logcounts
pval_ = unlist(lapply(rownames(logcounts_spe_cl30),function(xx){
  t.test(logcounts_spe_cl30[xx,which_m],logcounts_spe_cl30[xx,which_f])$p.value
}))
fc = unlist(lapply(rownames(logcounts_spe_cl30),function(xx){
  # t.test(logcounts_spe_cl30[xx,which_m],logcounts_spe_cl30[xx,which_f])$p.value
  mean(logcounts_spe_cl30[xx,which_m])-mean(logcounts_spe_cl30[xx,which_f])
}))

genes_deg = rownames(logcounts_spe_cl30)[pval_<.05 & abs(fc)>.2]
# sum(genes_deg)
# 255
genes_deg
######## plot, folder by cell cluster
library(ggpubr)
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

df_ = 
  # 


myvolcanoplot <- ggplot(data = df, aes(x = log2fc, y = -log10(pval), col = diffexpressed, label = delabel)) 
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') 
  geom_point(size = 2) 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable<br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)<br />
  coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits<br />
  labs(color = 'Severe', #legend_title,<br />
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis<br />
  ggtitle('Thf-like cells in severe COVID vs healthy patients') + # Plot title<br />
  geom_text_repel(max.overlaps = Inf) # To show all labels 



