# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaName="0"


###########
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)
aname <- "normcounts"


########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_nucleus.rds"))

se = spe_joint


########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female")
names(sex_perSampleID) = sampleID_unique

names(sex_perSampleID) = sampleID_unique

################################# find marker genes
######################################################################### 
# se <- scater::computeLibraryFactors(se)
# se <- scater::logNormCounts(se)
# sex_perSampleID
# se_seurat  <- as.Seurat(se)
# Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
# se_seurat = FindAllMarkers(se_seurat)
cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
# saveRDS(se_seurat, file = paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, "_nucleus.rds"))
sampleID_perCell = colData(se)[,"sample_id"]
sex = sex_perSampleID[as.character(sampleID_perCell)]
table(sex,sampleID_perCell)
table(sex,sampleID_perCell)

# se_seurat = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, "_nucleus.rds"))
######################################################################### 
### plot the results - heatmap for marker genes
######################################################################### 
# library(dplyr)
# genes_interest = c("gad1", "gad2", "slc32a1", "slc17a6",  "slc17a7")
# genes_interest = toupper(genes_interest)

######## process data into seurat format
obj = CreateSeuratObject(assays(se)$counts)
obj = NormalizeData(obj)
obj@meta.data$sex = sex
obj@meta.data$sampleID = sampleID_perCell
obj@meta.data$sex_sampleID = paste0(sex,"_",sampleID_perCell)

obj_clu15 = obj[,cellClu=="15"]
# obj_clu15 = FindVariableFeatures(obj_clu15)
VariableFeatures(obj_clu15) = rownames(procdata)
obj_clu15 = ScaleData(obj_clu15)
obj_clu15 = RunPCA(obj_clu15,npcs=5)
DimPlot(obj_clu15,reduction="pca",
        group.by = "sex")
pdf("featurePlot_clu15_ESR1_GAD1_SLC17A6_bySample.pdf", height = 30, width = 60)
FeaturePlot(obj_clu15,reduction="pca",
            features=c("GAD1","SLC17A6","ESR1"),
            split.by = "sex_sampleID",
            
            ncol=1)
dev.off()

table(sex[cellClu=="15"])
# Female   Male 
# 502   3630 

######## subset of cluster-15 cells
# counts_clu15 =  as.matrix(obj_clu15@assays$RNA@counts)
# dim(data_clu15)
# [1]  366 4132
# summary(counts_clu15["GAD1",])
# 0.00000 0.00000 0.00000 0.06412 0.00000 7.00000 
# 

######## 
counts_clu15 = assays(se)$counts[,cellClu == "15" ]
                                 # & colData(spe_joint)$sample_id=="colData(spe_joint)$sample_id"]

# counts_clu15 = counts_clu15[,1:1000]
exp_norm  = t(t(counts_clu15)/colSums(counts_clu15))*mean(colSums(counts_clu15))
exp_norm_kp  = exp_norm[which(rowVars(exp_norm)>0),]
# colSums(exp_norm_kp)


# > summary(counts_clu15["GAD1",])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   1.000   4.075   6.000  55.000 
# > # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   > # 0.000   0.000   1.000   4.075   6.000  55.000 
#   > summary(exp_norm["GAD1",])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   2.634   4.123   6.933  38.277 

######################################
######## pseudotime analysis of cluster 15 cells
######################################
library(igraph)
library(TSCAN)

procdata <- preprocess(data.matrix(exp_norm_kp),
                       minexpr_value = .2,
                       minexpr_percent = 0.2,
                       cvcutoff = .2)
# lpsmclust <- exprmclust(procdata,clusternum=3)#, clusternum=4
lpsmclust <- exprmclust(procdata)#, clusternum=4

plotmclust_v2(lpsmclust,
           show_cell_names=F,
           point_size = 0.02)



lpsorder
lpsorder <- TSCANorder(lpsmclust)

### plot genes 
GAD1expr <- log2(exp_norm_kp["GAD1",]+1)
singlegeneplot_v2(GAD1expr, TSCANorder(lpsmclust,
                                       flip=TRUE,
                                       orderonly=FALSE))
# 2767

SLC17A6expr <- log2(exp_norm_kp["SLC17A6",]+1)
singlegeneplot_v2(SLC17A6expr, TSCANorder(lpsmclust,
                                       flip=TRUE,
                                       orderonly=FALSE))
# 2767

expr <- log2(exp_norm_kp["ESR1",]+1)
singlegeneplot_v2(expr, TSCANorder(lpsmclust,
                                          flip=TRUE,
                                          orderonly=FALSE))

sex_clu15 = sex[cellClu=="15"]
singlegeneplot_v2(expr, TSCANorder(lpsmclust,
                                   flip=TRUE,
                                   orderonly=FALSE), subset_Cells =c(sex_clu15=="Male"),
                  ylim_= range(expr))

singlegeneplot_v2(expr, TSCANorder(lpsmclust,
                                   flip=TRUE,
                                   orderonly=FALSE), subset_Cells =c(sex_clu15=="Male"),
                  ylim_= range(expr))



singlegeneplot_v2(expr, TSCANorder(lpsmclust,
                                   flip=TRUE,
                                   orderonly=FALSE), subset_Cells =c(sex_clu15=="Female"),
                  ylim_= range(expr))
#################################
#################################
# 4191
sum(exp_norm_kp["GAD1",]>0)
# sum(GAD1expr>0)
# 323
# SLC17A6
sum(exp_norm_kp["SLC17A6",]>0)
sum(exp_norm_kp["SLC17A6",]>0&exp_norm_kp["GAD1",]>0)
# 29

summary(exp_norm_kp["GAD1",])

sum(exp_norm_kp["GAD1",]>0)

########################### 
# lpsorder <- TSCANorder(lpsmclust)
# lpsorder
plotmclust_v2 <- function (mclustobj, point_size = 0.05, x = 1, y = 2, MSTorder = NULL, show_tree = T, 
          show_cell_names = T, cell_name_size = 3, markerexpr = NULL) 
{
  color_by = "State"
  lib_info_with_pseudo <- data.frame(State = mclustobj$clusterid, 
                                     sample_name = names(mclustobj$clusterid))
  lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
  S_matrix <- mclustobj$pcareduceres
  pca_space_df <- data.frame(S_matrix[, c(x, y)])
  colnames(pca_space_df) <- c("pca_dim_1", "pca_dim_2")
  pca_space_df$sample_name <- row.names(pca_space_df)
  edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", 
                   by.y = "sample_name")
  edge_df$markerexpr <- markerexpr[edge_df$sample_name]
  if (!is.null(markerexpr)) {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2, 
                                    size = markerexpr))
    g <- g + geom_point(size = point_size,
                        aes_string(color = color_by), na.rm = TRUE)
  }
  else {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2))
    g <- g + geom_point( aes_string(color = color_by), na.rm = TRUE, 
                        size = point_size)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_tree) {
    clucenter <- mclustobj$clucenter[, c(x, y)]
    clulines <- NULL
    if (is.null(MSTorder)) {
      allsp <- shortest.paths(mclustobj$MSTtree)
      longestsp <- which(allsp == max(allsp), arr.ind = T)
      MSTorder <- get.shortest.paths(mclustobj$MSTtree, 
                                     longestsp[1, 1], longestsp[1, 2])$vpath[[1]]
    }
    for (i in 1:(length(MSTorder) - 1)) {
      clulines <- rbind(clulines, c(clucenter[MSTorder[i], 
      ], clucenter[MSTorder[i + 1], ]))
    }
    clulines <- data.frame(x = clulines[, 1], xend = clulines[, 
                                                              3], y = clulines[, 2], yend = clulines[, 4])
    g <- g + geom_segment(aes_string(x = "x", xend = "xend", 
                                     y = "y", yend = "yend", size = NULL), data = clulines, 
                          size = 1)
    clucenter <- data.frame(x = clucenter[, 1], y = clucenter[, 
                                                              2], id = 1:nrow(clucenter))
    g <- g + geom_text(aes_string(label = "id", x = "x", 
                                  y = "y", size = NULL), data = clucenter, size = 10)
  }
  g <- g + guides(colour = guide_legend(override.aes = list(size = 5))) + 
    xlab(paste0("PCA_dimension_", x)) + ylab(paste0("PCA_dimension_", 
                                                    y)) + theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(legend.position = "top", legend.key.size = unit(0.3, 
                                                          "in"), legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + 
    theme(legend.key = element_blank()) + theme(panel.background = element_rect(fill = "white")) + 
    theme(axis.text.x = element_text(size = 17, color = "darkred"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          axis.title.x = element_text(size = 20, vjust = -1), 
          axis.title.y = element_text(size = 20, vjust = 1), 
          plot.margin = unit(c(1, 1, 1, 1), "cm"))
  g
}

singlegeneplot_v2 <- function (geneexpr, TSCANorder, cell_size = 2,
                               subset_Cells=c(1:length(geneexpr)),
                               ylim_=NULL,
                               col_by=NULL) 
{
  Pseudotime <- NULL
  geneexpr <- geneexpr[TSCANorder[, 1]]
  exprdata <- cbind(TSCANorder, geneexpr)
  exprdata$State <- factor(exprdata$State)
  exprdata$predict <- fitted.values(mgcv::gam(geneexpr ~ s(Pseudotime, 
                                                           k = 3), data = exprdata))
  print(length(exprdata$Pseudotime))
  if(!is.null(col_by)){
    exprdata[,"group"] = col_by
    # q<-q+facet_wrap(~)
  }
  
  q <- ggplot(aes(Pseudotime, geneexpr), data = exprdata[subset_Cells,])
  q <- q + geom_point(aes_string(color = "State"), size = I(cell_size))
  q <- q + geom_line(aes(Pseudotime, predict), data = exprdata)
  q <- q + ylab("Expression") + xlab("Pseudotime")
  q <- q + theme(strip.background = element_rect(colour = "white", 
                                                 fill = "white")) + 
    theme(panel.border = element_blank(), 
          axis.line = element_line()) + theme(panel.grid.minor.x = element_blank(), 
                                              panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.major.y = element_blank()) + theme(panel.background = element_rect(fill = "white"))

  if(!is.null(ylim_)){
          q <- q + ylim(ylim_)
  }
  if(!is.null(col_by)){
    
    q<-q+facet_wrap(~group,ncol=2)
  }
  q
}


