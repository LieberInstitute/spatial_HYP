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
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)
dir_data_processed =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"

########################################################################
########### load ARC annotation result
knn_pred_arc_c_ = readRDS("knn_ARC_twoRounds.rds")
knn_pred_vmh_c_ = readRDS("knn_VMH_twoRounds.rds")

###################################################################### 
################################### load multi-sample Banksy result
###################################################################### 
############ load cell type clustering result, nucleus only
lambda <- 0
lambdaName="0"
se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", 2, ".rds"))

cellClu_annotation = colData(se)[,"clust_M0_lam0_k50_res2"]
length(cellClu_annotation)

########################################################################
################################# format data
######################################################################### 
se <- scater::computeLibraryFactors(se)
##
# which_cells_kp =  which( knn_pred_arc_c_ =="ARC")

##
cellClu_domain_kp = knn_pred_arc_c_
cellClu_annotation_kp =as.character(cellClu_annotation)
##
se <- scater::logNormCounts(se[,])

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]


######## load DEG list for each ARC cell type

list_DEG_perClu = readRDS("list_DEG_perClu.rds")
# saveRDS(list_DEG_perClu, file = "list_DEG_perClu.rds")
# genes_deg

### assign levels to the factor denoting sample ID based on sex - male and female (for better organized in the plotting)

################################################################################ 
################################ plot, folder by cell cluster
################################################################################ 
library(ggpubr)
gene_tmp = "GLRA2"
  
df <- as.data.frame(colData(se))
        # ( logcounts(se)[gene_tmp, ])))
df$expr = logcounts(se)[gene_tmp, ]

kp_ = df$sample_id!="br8667c"
df = df[kp_,]

SampleIDs_kp = unique(df$sample_id)
nsample_kp = length(SampleIDs_kp)
df$sample_id = factor(df$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
df$sex = factor(df$sex, level = c("Male", "Female"))


# ARC
pseudobulk_clu_tmp =  unlist(lapply(SampleIDs_kp,function(xx){
  mean(df[which(df$sample_id == xx & knn_pred_arc_c_[kp_]=="ARC"), "expr"])
}))

pseudobulk_clu_tmp_back =  unlist(lapply(SampleIDs_kp,function(xx){
  mean(df[which(df$sample_id == xx & knn_pred_arc_c_[kp_]!="ARC"), "expr"])
}))

df_pseudobulk_tmp_c = data.frame(expr = c(pseudobulk_clu_tmp, pseudobulk_clu_tmp_back),
                               sample_id = rep(SampleIDs_kp, 2),
                               sex = sex_perSampleID[rep(SampleIDs_kp, 2)],
                               group = rep(c("within ARC","outside ARC"), each = length(SampleIDs_kp)))

df_pseudobulk_tmp_c$sample_id = factor(df_pseudobulk_tmp_c$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
df_pseudobulk_tmp_c$sex = factor(df_pseudobulk_tmp_c$sex, level = c("Male", "Female"))
df_pseudobulk_tmp_c$group = factor(df_pseudobulk_tmp_c$group, level = c("within ARC","outside ARC"))

pdf("pseudobulk_ValidateARC.pdf")
p12_pseudobulk <- ggplot(df_pseudobulk_tmp_c, aes(x = group, y = expr, fill = group)) + 
  theme_bw() +
  geom_violin(width=1,position=position_dodge(1)) +
  geom_boxplot(width=0.1, position=position_dodge(1)) +
  geom_point() +
  theme(plot.title = element_text(face = "italic"),
        panel.grid = element_blank()
  ) + facet_wrap(~ sex, ncol = 2) +
  scale_fill_manual(values=c( "#56B4E9","#E69F00" ))
p12_pseudobulk
dev.off()

# VMH
gene_tmp = "ANKRD34B"
df$expr  = logcounts(se)[gene_tmp, ]

pseudobulk_clu_tmp =  unlist(lapply(SampleIDs_kp,function(xx){
  mean(df[which(df$sample_id == xx & knn_pred_vmh_c_[kp_]=="VMH"), "expr"])
}))

pseudobulk_clu_tmp_back =  unlist(lapply(SampleIDs_kp,function(xx){
  mean(df[which(df$sample_id == xx & knn_pred_vmh_c_[kp_]!="VMH"), "expr"])
}))

df_pseudobulk_tmp_c = data.frame(expr = c(pseudobulk_clu_tmp, pseudobulk_clu_tmp_back),
                                 sample_id = rep(SampleIDs_kp, 2),
                                 sex = sex_perSampleID[rep(SampleIDs_kp, 2)],
                                 group = rep(c("within VMH","outside VMH"), each = length(SampleIDs_kp)))

df_pseudobulk_tmp_c$sample_id = factor(df_pseudobulk_tmp_c$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
df_pseudobulk_tmp_c$sex = factor(df_pseudobulk_tmp_c$sex, level = c("Male", "Female"))
df_pseudobulk_tmp_c$group = factor(df_pseudobulk_tmp_c$group, level = c("within VMH","outside VMH"))

pdf("pseudobulk_ValidateVMH.pdf")
p12_pseudobulk <- ggplot(df_pseudobulk_tmp_c, aes(x = group, y = expr, fill = group)) + 
  theme_bw() +
  geom_violin(width=1,position=position_dodge(1)) +
  geom_boxplot(width=0.1, position=position_dodge(1)) +
  geom_point() +
  theme(plot.title = element_text(face = "italic"),
        panel.grid = element_blank()
  ) + facet_wrap(~ sex, ncol = 2) +
  scale_fill_manual(values=c( "#56B4E9","#E69F00" ))
p12_pseudobulk
dev.off()







