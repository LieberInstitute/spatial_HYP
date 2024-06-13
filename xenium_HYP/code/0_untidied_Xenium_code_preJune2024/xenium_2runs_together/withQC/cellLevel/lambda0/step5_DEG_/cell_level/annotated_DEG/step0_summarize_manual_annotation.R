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

dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)

############ load cell type clustering result, cell level
lambda <- 0
lambdaName="0"
se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", 2, ".rds"))



####### ARC DEG
cellClu_annotation = colData(se)[,"clust_M0_lam0_k50_res2"]

nclu= length(unique(cellClu_annotation))
####### put manual annotation
cellTypeAnnotations = as.character(1:nclu)
names(cellTypeAnnotations) = as.character(1:nclu)
cellTypeAnnotations[1] = "Other_GABAergic"
cellTypeAnnotations[2] = "Oligodendrocyte_Precursor_cells_(OPC)"
cellTypeAnnotations[3] = "Insufficient_markers"
cellTypeAnnotations[4] = "Astrocyte"
cellTypeAnnotations[5] = "Endothelium"
cellTypeAnnotations[6] = "Oligodendrocyte"
cellTypeAnnotations[7] = "Astrocyte_2"
cellTypeAnnotations[8] = "Oligodendrocyte"
cellTypeAnnotations[9] = "Vascular_uncertain"
cellTypeAnnotations[10] = "Oligodendrocyte_2_(OPALIN+)"
cellTypeAnnotations[11] = "Oligodendrocyte_3_(OPALIN+)"
cellTypeAnnotations[12] = "ARC_1_(POMC)"
cellTypeAnnotations[13] = "Microglia_1"
cellTypeAnnotations[14] = "Excitatory_Periventricular_(CRH,_TRH)"
cellTypeAnnotations[15] = "Tanycyte_1"
cellTypeAnnotations[16] = "Microglia_(Reactive)"
cellTypeAnnotations[17] = "Mixed_Astrocyte-Arc_POMC"
cellTypeAnnotations[18] = "VMH_1_(Excitatory)"
cellTypeAnnotations[19] = "VMH_2_(Excitatory)"
cellTypeAnnotations[20] = "Unsure1"
cellTypeAnnotations[21] = "DISCARD"
cellTypeAnnotations[22] = "Supraoptic_nu_(AVP+OXT+)_1"
cellTypeAnnotations[23] = "VMH_3"
cellTypeAnnotations[24] = "Macrophage"
cellTypeAnnotations[25] = "ARC_2_(GHRH_GAL)"
cellTypeAnnotations[26] = "Oligodendrocyte_4_(Optic_Nerve)"
cellTypeAnnotations[27] = "Microglia_2"
cellTypeAnnotations[28] = "Tanycyte_2"
cellTypeAnnotations[29] = "Peripheral_immune_NK-or_Tcell"
cellTypeAnnotations[30] = "Astrocyte_3"
cellTypeAnnotations[31] = "Supraoptic_nu_(AVP+OXT+)_2"
cellTypeAnnotations[32] = "VMH_Lateral_border"
cellTypeAnnotations[33] = "ARC_3_(TAC3_ESR1)"
cellTypeAnnotations[34] = "VMH_4"
cellTypeAnnotations[35] = "DISCARD"
cellTypeAnnotations[36] = "VMH_5"
cellTypeAnnotations[37] = "DISCARD"
cellTypeAnnotations[38] = "DISCARD"
cellTypeAnnotations[39] = "Mixed_Oligodendrocytes_and_Supraoptic_nu_(AVP+OXT+)_3"
cellTypeAnnotations[40] = "VMH_6_Lateral_bit"
cellTypeAnnotations[41] = "GABAergic_unspecified"

saveRDS(cellTypeAnnotations, file = "cellTypeAnnotations.rds")
