# ml conda_R/4.3
# R
###########
# library(Seurat)
# library(SummarizedExperiment)
# library(SpatialExperiment)
# library(Banksy)
library(scuttle)
# library(scater)
library(cowplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)

###########
list_pval_allGenes_perClu = readRDS("list_pval_allGenes_perClu_VMH.rds")
list_fc_allGenes_perClu = readRDS("list_fc_allGenes_perClu_VMH.rds")
list_DEG_perClu = readRDS("list_DEG_perClu_VMH.rds")
# genes_deg_tmp = c("LYPD6B", "CDH13","MC3R","CDH4","RXRG")
# genes_deg_tmp = c("LYPD6B", "CDH13","CDH4")

mat_pval_allGenes_perClu = matrix(unlist(list_pval_allGenes_perClu), ncol = length(list_pval_allGenes_perClu))
colnames(mat_pval_allGenes_perClu) = names(list_pval_allGenes_perClu)
rownames(mat_pval_allGenes_perClu) = names(list_pval_allGenes_perClu[[1]])


mat_FC_allGenes_perClu = matrix(unlist(list_fc_allGenes_perClu), ncol = length(list_fc_allGenes_perClu))
colnames(mat_FC_allGenes_perClu) = names(list_fc_allGenes_perClu)
rownames(mat_FC_allGenes_perClu) = names(list_fc_allGenes_perClu[[1]])

mat_pval_allGenes_perClu = mat_pval_allGenes_perClu + 10^ (-10)


rowMins_pval = rowMins(abs(mat_pval_allGenes_perClu))
rowMaxs_fc = rowMaxs(abs(mat_FC_allGenes_perClu))

kp_genes_which = which(rowMins_pval<0.05 & rowMaxs_fc> 0.5)

length(kp_genes_which)

########### plot p-value
library(pheatmap)

mat_pval_unsigned =  -log10(mat_pval_allGenes_perClu[kp_genes_which, ] )
mat_pval_signed = mat_pval_unsigned
mat_pval_signed[mat_FC_allGenes_perClu[kp_genes_which,] < 0] = -mat_pval_signed[mat_FC_allGenes_perClu[kp_genes_which,] < 0] 

set.seed(111)
clu_genes = kmeans(mat_pval_signed, 2)$cluster
set.seed(111)
clu_cellClus = kmeans(t(mat_pval_signed), 2)$cluster

# mat_pval_signed[mat_FC_allGenes_perClu[kp_genes_which,] < 0] = - mat_pval_signed[mat_FC_allGenes_perClu[kp_genes_which,] < 0]
mat_pval_signed_MaleHigher = mat_pval_unsigned
mat_pval_signed_MaleHigher[mat_FC_allGenes_perClu[kp_genes_which,] < 0] = min(mat_pval_unsigned)
mat_pval_signed_FemaleHigher = mat_pval_unsigned
mat_pval_signed_FemaleHigher[mat_FC_allGenes_perClu[kp_genes_which,] > 0] = min(mat_pval_unsigned)

pdf(paste0("representative_plots/heatmap_pval_perCellClu_signed_VMH.pdf"),width=15,height=15)
p  = pheatmap(mat_pval_signed_MaleHigher[order(clu_genes),order(clu_cellClus)], cluster_rows = F, cluster_cols = F) #, cluster_rows = F, cluster_cols = F, annotation_col = col_annotation, annotation_row= row_annotation,fontsize_row=1
print(p)
p  = pheatmap(mat_pval_signed_FemaleHigher[order(clu_genes),order(clu_cellClus)], cluster_rows = F, cluster_cols = F) #, cluster_rows = F, cluster_cols = F, annotation_col = col_annotation, annotation_row= row_annotation,fontsize_row=1
print(p)
dev.off()


########### plot log fc
mat_FC_allGenes_perClu_kp = mat_FC_allGenes_perClu[kp_genes_which,]
mat_FC_allGenes_perClu_kp[mat_pval_allGenes_perClu[kp_genes_which, ]>0.05] = NA


pdf(paste0("representative_plots/heatmap_logFC_perCellClu_signed_VMH.pdf"),width=15,height=15)
p  = pheatmap(mat_FC_allGenes_perClu_kp[order(clu_genes),order(clu_cellClus)], cluster_rows = F, cluster_cols = F) #, cluster_rows = F, cluster_cols = F, annotation_col = col_annotation, annotation_row= row_annotation,fontsize_row=1
print(p)
# p  = pheatmap(mat_pval_signed_FemaleHigher[order(clu_genes),order(clu_cellClus)], cluster_rows = F, cluster_cols = F) #, cluster_rows = F, cluster_cols = F, annotation_col = col_annotation, annotation_row= row_annotation,fontsize_row=1
# print(p)
dev.off()





