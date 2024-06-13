# ml conda_R/4.3
# R

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
spe_joint=se
if_ARC = array("not ARC", dim=ncol(se))
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

########### keep cells that are in ARC domain and are kept in both process - full data and nucleo
# cellIDs_kp = intersect(cellID_cellClustAnnotation,
                       # cellIDs_se_domain[cellClu_domain=="3"])
# length(cellIDs_kp)
# 366887
# 879396
########### 

# cellClu_domain_kp = colData(se_domain)[cellIDs_kp,paste0("clust_M0_lam", 1, "_k50_res0.6")]
# cellClu_annotation_kp = colData(se)[cellIDs_kp,paste0("clust_M0_lam", 0, "_k50_res0.6")]
# table(cellClu_domain_kp,cellClu_annotation_kp)


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
####### ARC DEG
library(ggpubr)
library(ggplot2)
library(cowplot)


spe = se[,]
# cellType_clu_nuc_used = cellClu_annotation_kp
cellType_clu_unique= sort(unique(cellClu_annotation_kp))
# > dim(genes_ARC_sexDEG)
# [1] 15  7
# > 
cellTypeAnnotations = cellClu_annotation_kp
# cellTypeAnnotations = c("astrocytes", "likely more than one cell types","myelinating oligodendrocytes",
#                         "microglia","","oligodendrocyte precursor",
#                         "immune related","epithelial cells","astrocytes",
#                         "excitatory neurons","", "POMC neuronal cells",
#                         "non-myelinating matureish(?) oligodendrocytes","epithelial cells","",
#                         "","T cells","AVP neurons",
#                         "monocyte")
names(cellTypeAnnotations) = as.character(1:length(cellType_clu_unique))

pdf("AnnotatedSexDEG_ARC_boxplot_perCellCluster.pdf",width=30,height=30)
for(gene_tmp in genes_ARC_sexDEG[,1] ){
  
  plist = list()
  
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  df = df[df$sample_id!="br8667c",]
  
  
  for(cellType_clu in cellType_clu_unique){
    
    
    df_cellType_clu_tmp = df[cellClu_annotation_kp == cellType_clu,]
    df_cellType_clu_tmp$sample_id = factor(df_cellType_clu_tmp$sample_id, level = sample_names[order(sex_perSampleID)])
    # sample_names_new = c("br1225a","br1225b","br8741c",
                         # "br8741d","br5459a","br5459b",
                         # "br8667c")
    # sample_names= c(sample_names_old, sample_names_new)
    # sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female",
                        # c("Male","Male","Female","Female","Male","Male","Female"))
    plist[[as.integer(cellType_clu)]]  <- ggplot(df_cellType_clu_tmp, aes(x = sex, y = expr, fill = sample_id)) + 
      ggtitle(paste0(gene_tmp,", cluster ", cellType_clu)) + 
      theme_bw() + 
      geom_violin(width=1, position = position_dodge(1)) +
      geom_boxplot(width=0.1, position = position_dodge(1)) +
      theme(plot.title = element_text(face = "italic"),
            panel.grid = element_blank()
      )+
      scale_fill_manual(values=c( "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#E69F00","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00" ))
    
  }

  
  print(plot_grid(plotlist=plist, ncol=4))
  # print(plot_grid(plist[[1]],plist[[2]], ncol=4))
  
}
dev.off()



####################################### summarize into a gene by cell-type table
ncluster = length(cellType_clu_unique)
mat_pvalue_F_higher_than_M = mat_pvalue_F_lower_than_M=array(dim=c(ncluster,15))

colnames(mat_pvalue_F_higher_than_M) = colnames(mat_pvalue_F_lower_than_M) = genes_ARC_sexDEG[,1]
rownames(mat_pvalue_F_higher_than_M) = rownames(mat_pvalue_F_lower_than_M) = paste0("cluster ",1:ncluster)

for(gene_tmp in genes_ARC_sexDEG[,1] ){
  
  df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
  df$sample_id=colData(spe)$sample_id
  
  df$sex=sex_perSampleID[as.character(df$sample_id)]
  df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
  
  for(cellType_clu in cellType_clu_unique){
    df_cellType_clu_tmp = df[cellClu_annotation_kp == cellType_clu,]
    # plist[[as.integer(cellType_clu)]]  <- ggplot(df_cellType_clu_tmp, aes(x = sex, y = expr, fill = sample_id)) + 
    which_m = which(df_cellType_clu_tmp$sex=="Male")
    which_f = which(df_cellType_clu_tmp$sex=="Female")
    
    mat_pvalue_F_higher_than_M[paste0("cluster ",cellType_clu), gene_tmp] = t.test(df_cellType_clu_tmp$expr[which_m],
                                                                                   df_cellType_clu_tmp$expr[which_f],
                                                                                   alternative = "less")$p.value
    
    
    mat_pvalue_F_lower_than_M[paste0("cluster ",cellType_clu), gene_tmp] = t.test(df_cellType_clu_tmp$expr[which_m],
                                                                                   df_cellType_clu_tmp$expr[which_f],
                                                                                   alternative = "greater")$p.value
  }
}

save(mat_pvalue_F_higher_than_M,mat_pvalue_F_lower_than_M,file="pval_DEG_ARC_perCellType.RData")

rownames(mat_pvalue_F_higher_than_M) = rownames(mat_pvalue_F_lower_than_M) = paste0(paste0("cluster ",1:ncluster))
write.csv(round(mat_pvalue_F_higher_than_M,2), file="mat_pvalue_F_higher_than_M_perCellType_ARC.csv",quote=F)
write.csv(round(mat_pvalue_F_lower_than_M,2), file="mat_pvalue_F_lower_than_M_perCellType_ARC.csv",quote=F)

mat_pvalue_F_higher_than_M_bin =mat_pvalue_F_lower_than_M_bin = array("",dim=dim(mat_pvalue_F_higher_than_M))
mat_pvalue_F_higher_than_M_bin[mat_pvalue_F_higher_than_M<=0.05] = "X"
mat_pvalue_F_lower_than_M_bin[mat_pvalue_F_lower_than_M<=0.05] = "X"

dimnames(mat_pvalue_F_higher_than_M_bin) = dimnames(mat_pvalue_F_lower_than_M_bin) = dimnames(mat_pvalue_F_higher_than_M)
write.csv(mat_pvalue_F_higher_than_M_bin, file="mat_pvalue_F_higher_than_M_perCellType_ARC_binary.csv",quote=F)
write.csv(mat_pvalue_F_lower_than_M_bin, file="mat_pvalue_F_lower_than_M_perCellType_ARC_bin.csv",quote=F)

mat_bin_overall= array("",dim=dim(mat_pvalue_F_higher_than_M))
dimnames(mat_bin_overall)=dimnames(mat_pvalue_F_higher_than_M_bin) 
mat_bin_overall[mat_pvalue_F_higher_than_M_bin=="X"]="F > M"
mat_bin_overall[mat_pvalue_F_lower_than_M_bin=="X"]="F < M"

write.csv(mat_bin_overall, file="mat_mat_bin_overall_perCellType_ARC_overall.csv",quote=F)

########################################################## 
################ adjust the p-value
########################################################## 
# dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# setwd(dir_out)

load("pval_DEG_ARC_perCellType.RData")
# save(mat_pvalue_F_higher_than_M,mat_pvalue_F_lower_than_M,file="pval_DEG_ARC_perCellType.RData")
mat_pvalue_F_lower_than_M_Adj = mat_pvalue_F_higher_than_M_Adj = array(dim=dim(mat_pvalue_F_higher_than_M))
nTest = length(mat_pvalue_F_lower_than_M_Adj)

mat_pvalue_F_higher_than_M_Adj[1:nTest] = p.adjust(mat_pvalue_F_higher_than_M[1:nTest] , 
                                                   method = "BH", n = nTest)
mat_pvalue_F_lower_than_M_Adj[1:nTest] =  p.adjust(mat_pvalue_F_lower_than_M[1:nTest], 
                                                   method = "BH", n = nTest)

save(mat_pvalue_F_higher_than_M_Adj,mat_pvalue_F_lower_than_M_Adj,
     file="pval_DEG_ARC_perCellTyp_Adj.RData")


dimnames(mat_pvalue_F_higher_than_M_Adj) = dimnames(mat_pvalue_F_lower_than_M_Adj) = dimnames(mat_pvalue_F_higher_than_M)
write.csv(round(mat_pvalue_F_higher_than_M_Adj,2), file="mat_pvalue_F_higher_than_M_Adj_perCellType_ARC.csv",quote=F)
write.csv(round(mat_pvalue_F_lower_than_M_Adj,2), file="mat_pvalue_F_lower_than_M_Adj_perCellType_ARC.csv",quote=F)

mat_pvalue_F_higher_than_M_adj_bin =mat_pvalue_F_lower_than_M_adj_bin = array("",dim=dim(mat_pvalue_F_higher_than_M))
mat_pvalue_F_higher_than_M_adj_bin[mat_pvalue_F_higher_than_M_Adj<=0.05] = "X"
mat_pvalue_F_lower_than_M_adj_bin[mat_pvalue_F_lower_than_M_Adj<=0.05] = "X"


dimnames(mat_pvalue_F_higher_than_M_adj_bin) = dimnames(mat_pvalue_F_lower_than_M_adj_bin) = dimnames(mat_pvalue_F_higher_than_M)
write.csv(mat_pvalue_F_higher_than_M_adj_bin, file="mat_pvalue_F_higher_than_M_perCellType_ARC_adj_binary.csv",quote=F)
write.csv(mat_pvalue_F_lower_than_M_adj_bin, file="mat_pvalue_F_lower_than_M_perCellType_ARC_adj_bin.csv",quote=F)



mat_adj_bin_overall= array("",dim=dim(mat_pvalue_F_higher_than_M_Adj))
dimnames(mat_adj_bin_overall)=dimnames(mat_pvalue_F_higher_than_M_adj_bin) 
mat_adj_bin_overall[mat_pvalue_F_higher_than_M_adj_bin=="X"]="F > M"
mat_adj_bin_overall[mat_pvalue_F_lower_than_M_adj_bin=="X"]="F < M"

write.csv(mat_adj_bin_overall, file="mat_adj_bin_perCellType_ARC_overall.csv",quote=F)


################  add FC to the table
mat_adj_bin_overall = read.csv("mat_adj_bin_perCellType_ARC_overall.csv")


load("FC_DEG_ARC_perCellType.RData")
# save(mat_FC_F_higher_than_M,mat_FC_F_lower_than_M,file="FC_DEG_ARC_perCellType.RData")

mat_adj_bin_overall_withFC = as.matrix(mat_adj_bin_overall[,-1])
which_F_higher = which(mat_adj_bin_overall_withFC=="F > M" & mat_FC_F_higher_than_M>1.5)
which_M_higher = which(mat_adj_bin_overall_withFC=="F < M" & mat_FC_F_lower_than_M>1.5)

# mat_log2FC_significant = mat_adj_bin_overall_withFC
# dimnames(mat_log2FC_significant) = dimnames((mat_adj_bin_overall_withFC))
mat_adj_bin_overall_withFC[which_F_higher] = paste0(mat_adj_bin_overall_withFC[which_F_higher], 
                                                    " (FC=",
                                                    round(mat_FC_F_higher_than_M[which_F_higher],1),
                                                    ")")
# mat_adj_bin_overall_withFC[which_M_higher] = mat_log2FC_F_lower_than_M[which_M_higher]
mat_adj_bin_overall_withFC[which_M_higher] =paste0(mat_adj_bin_overall_withFC[which_M_higher], 
                                                   " (FC=",
                                                   round(mat_FC_F_lower_than_M[which_M_higher],1),
                                                   ")")


mat_adj_bin_overall_withFC[-c(which_M_higher,which_F_higher )] = ""

# rownames(mat_adj_bin_overall_withFC) = rownames(mat_FC_F_lower_than_M)

# write.csv(mat_adj_bin_overall_withFC, file="mat_adj_bin_overall_withFC_ARC.csv",quote=F)


# cellTypeAnnotations = readRDS("cellTypeAnnotations.rds")
# rownames(mat_adj_bin_overall_withFC) = cellTypeAnnotations

# rownames(mat_adj_bin_overall_withFC) = paste0(rownames(mat_FC_F_lower_than_M) ,"----", cellTypeAnnotations) 

write.csv(mat_adj_bin_overall_withFC, file="mat_adj_bin_overall_withFC_ARC.csv",quote=F)

# VMH: VSNL
# ARC: CRH,SSTR