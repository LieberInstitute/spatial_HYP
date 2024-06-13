
##########
lambdaName=0
lambda=0
# res <- 2
###########
# dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
# dir_data_processed =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
# dir_data_processed_banksy =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# setwd(dir_out)
###########
# dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun_nucleus/"
# dir.create(dir_data_processed_banksy)
# dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out_nucleus/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")
setwd(dir_out)

########### load cell type clustering using nucle onlu
# cellType_clu_nuc = readRDS("cellType_clu_nuc.rds")
# saveRDS(cellType_clu_nuc, file = "cellType_clu_nuc.rds")
###########

# geneInfo_ = read.csv(paste0(dir_data_processed,"data/Xenium_HYP_Candidate_Genes.csv"))
# head(geneInfo_)

###########
########### find genes that are annotated
###########

# geneInfo_[,2]

# genes_ARC = geneInfo_[grep("ARC",geneInfo_[,2]),]
# genes_ARC_sexDEG = genes_ARC[grep("Sex-differential",genes_ARC[,2]),]
dir_data_processed =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"

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
########### find genes that are annotated
###########
geneInfo_ = read.csv(paste0(dir_data_processed,"/Xenium_HYP_Candidate_Genes.csv"))
head(geneInfo_)

geneInfo_[,2]

genes_ARC = geneInfo_[grep("ARC",geneInfo_[,2]),]
genes_ARC_marker1 = genes_ARC[grep("ARC Enriched",genes_ARC[,2]),]
genes_ARC_marker2 = genes_ARC[grep("Known ARC marker",genes_ARC[,2]),]
genes_ARC_marker3 = genes_ARC[grep("ARC enriched",genes_ARC[,2]),]
genes_ARC_marker = rbind(genes_ARC_marker1, genes_ARC_marker2, genes_ARC_marker3)
genes_ARC_marker[,2]


genes_ARC = geneInfo_[grep("ARC",geneInfo_[,2]),]
genes_ARC_sexDEG = genes_ARC[grep("Sex-differential",genes_ARC[,2]),]
genes_ARC_sexDEG[,2]

########### load multi-sample Banksy result
############ load domain segmentation result
lambda <- 1
lambdaName="1"
# spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
# cellClu_domain = colD

############ load cell type clustering result, nucleus only
lambda <- 0
lambdaName="0"
se = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", 2, ".rds"))
# se = spe_joint
# cellID_cellClustAnnotation=colData(se)$cell_id

cellClu_annotation = colData(se)[,"clust_M0_lam0_k50_res2"]
length(cellClu_annotation)
# 879396
########### get ARC domain
if_ARC = array("not ARC", dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("/dcs04/hansen/data/ywang/xenium_joint_runs/out/knn_predARC_",sample,"_k500_2ndRun.rds"))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
  if_ARC[which_sample][knn_pred>0.3] = "ARC"
}
# colData(spe_joint)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
# se_domain = spe_joint
# cellIDs_se_domain = colData(se_domain)$cell_id
length(if_ARC)
# 879396


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
################################# format data
######################################################################### 

se <- scater::computeLibraryFactors(se)
##
which_cells_kp =  which(colData(se)$sizeFactor>0)

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
library(ggpubr)
library(ggplot2)
library(cowplot)



spe = se[,]
# cellType_clu_nuc_used = cellClu_annotation_kp
cellType_clu_unique= sort(unique(cellClu_annotation_kp))
dim(genes_ARC_sexDEG[-1,])
# cellTypeAnnotations = cellClu_annotation_kp
# names(cellTypeAnnotations) = as.character(1:length(cellType_clu_unique))

# [1] 15  7

genes_ARC_sexDEG[-1,1]%in%rownames(spe)
# [1]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
# [13]  TRUE  TRUE  TRUE

genes_ARC_sexDEG_kp = intersect(genes_ARC_sexDEG[-1,1],
                                rownames(spe))

# cellTypeAnnotations = c("astrocytes", "likely more than one cell types","myelinating oligodendrocytes",
#                         "microglia","","oligodendrocyte precursor",
#                         "immune related","epithelial cells","astrocytes",
#                         "excitatory neurons","", "POMC neuronal cells",
#                         "non-myelinating matureish(?) oligodendrocytes","epithelial cells","",
#                         "","T cells","AVP neurons",
#                         "monocyte")
names(cellTypeAnnotations) = as.character(1:ncluster)
names(cellTypeAnnotations) = as.character(1:length(cellType_clu_unique))

# saveRDS(cellTypeAnnotations, file = "cellTypeAnnotations.rds")
# write.csv(cellTypeAnnotations, file="cellTypeAnnotations.csv",quote=F)

####################################### summarize into a gene by cell-type table
mat_log2FC_F_higher_than_M = mat_log2FC_F_lower_than_M=array(dim=c(ncluster,14))

colnames(mat_log2FC_F_higher_than_M) = colnames(mat_log2FC_F_lower_than_M) = genes_ARC_sexDEG_kp
rownames(mat_log2FC_F_higher_than_M) = rownames(mat_log2FC_F_lower_than_M) = paste0("cluster ",1:ncluster)

for(gene_tmp in genes_ARC_sexDEG_kp ){
  # counts(spe)
  df <- as.data.frame(
    cbind(spatialCoords(spe),
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id

  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)

  for(cellType_clu in cellType_clu_unique){
    df_cellType_clu_tmp = df[cellClu_annotation_kp == cellType_clu,]

    which_m = which(df_cellType_clu_tmp$sex=="M")
    which_f = which(df_cellType_clu_tmp$sex=="F")

    mat_log2FC_F_higher_than_M[paste0("cluster ",cellType_clu), gene_tmp] = mean(df_cellType_clu_tmp$expr[which_f]) - mean(df_cellType_clu_tmp$expr[which_m])
                                                                                  # t.test(df_cellType_clu_tmp$expr[which_m],
                                                                                  # df_cellType_clu_tmp$expr[which_f],
                                                                                  # alternative = "less")$p.value

    mat_log2FC_F_lower_than_M[paste0("cluster ",cellType_clu), gene_tmp] = mean(-df_cellType_clu_tmp$expr[which_f]) + mean(df_cellType_clu_tmp$expr[which_m])
  }
}


# 
save(mat_log2FC_F_higher_than_M,mat_log2FC_F_lower_than_M,file="log2FC_DEG_ARC_perCellType.RData")
# 
rownames(mat_log2FC_F_higher_than_M) = rownames(mat_log2FC_F_lower_than_M) = paste0(paste0("cluster ",1:ncluster))
write.csv(round(mat_log2FC_F_higher_than_M,2), file="mat_log2FC_F_higher_than_M_perCellType_ARC.csv",quote=F)
write.csv(round(mat_log2FC_F_lower_than_M,2), file="mat_log2FC_F_lower_than_M_perCellType_ARC.csv",quote=F)
genes_ARC_sexDEG_kp =genes_ARC_sexDEG[,1]
mat_FC_F_higher_than_M = mat_FC_F_lower_than_M=array(dim=c(ncluster,15))

colnames(mat_FC_F_higher_than_M) = colnames(mat_FC_F_lower_than_M) = genes_ARC_sexDEG_kp
rownames(mat_FC_F_higher_than_M) = rownames(mat_FC_F_lower_than_M) = paste0("cluster ",1:ncluster)

for(gene_tmp in genes_ARC_sexDEG[,1]  ){
  # counts(spe)
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = counts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  
  for(cellType_clu in cellType_clu_unique){
    df_cellType_clu_tmp = df[cellClu_annotation_kp == cellType_clu,]
    
    which_m = which(df_cellType_clu_tmp$sex=="Male")
    which_f = which(df_cellType_clu_tmp$sex=="Female")
    
    mat_FC_F_higher_than_M[paste0("cluster ",cellType_clu), gene_tmp] = mean(df_cellType_clu_tmp$expr[which_f])/mean(df_cellType_clu_tmp$expr[which_m])
    # t.test(df_cellType_clu_tmp$expr[which_m],
    # df_cellType_clu_tmp$expr[which_f],
    # alternative = "less")$p.value
    
    mat_FC_F_lower_than_M[paste0("cluster ",cellType_clu), gene_tmp] =    mean(df_cellType_clu_tmp$expr[which_m]) / mean(df_cellType_clu_tmp$expr[which_f])
  }
}

save(mat_FC_F_higher_than_M,mat_FC_F_lower_than_M,file="FC_DEG_ARC_perCellType.RData")

rownames(mat_FC_F_higher_than_M) = rownames(mat_FC_F_lower_than_M) = paste0(paste0("cluster ",1:ncluster))
write.csv(round(mat_FC_F_higher_than_M,2), file="mat_FC_F_higher_than_M_perCellType_ARC.csv",quote=F)
write.csv(round(mat_FC_F_lower_than_M,2), file="mat_FC_F_lower_than_M_perCellType_ARC.csv",quote=F)

# ####################################### 


