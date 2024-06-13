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
which_cells_kp =  which(colData(se)$sizeFactor>0 & knn_pred_arc_c_ =="ARC")

##
cellClu_domain_kp = knn_pred_arc_c_[which_cells_kp]
cellClu_annotation_kp =as.character(cellClu_annotation[which_cells_kp])
##
se <- scater::logNormCounts(se[,which_cells_kp])

se_seurat  <- as.Seurat(se)
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]


######## load DEG list for each ARC cell type

list_DEG_perClu = readRDS("list_DEG_perClu.rds")
# saveRDS(list_DEG_perClu, file = "list_DEG_perClu.rds")
# genes_deg

################################################################################ 
################################ plot, folder by cell cluster
################################################################################ 
library(ggpubr)

df <- as.data.frame(
  cbind(spatialCoords(se), 
       t( logcounts(se)[unique(unlist(list_DEG_perClu)), ])))


df$sample_id=colData(se)$sample_id

df$sex=sex_perSampleID[as.character(df$sample_id)]
df$sex_sampleID=paste0(df$sex,"_",df$sample_id)

### remove one sample where it seems has little ARC contained
df = df[df$sample_id!="br8667c",]

### assign levels to the factor denoting sample ID based on sex - male and female (for better organized in the plotting)
SampleIDs_kp = unique(df$sample_id)
nsample_kp = length(SampleIDs_kp)
df$sample_id = factor(df$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
df$sex = factor(df$sex, level = c("Male", "Female"))
table(df$sex)


#####
dir.create("representative_plots")
genes_deg_tmp = c("LYPD6B", "CDH13","MC3R","CDH4","RXRG")
for(gene_tmp in genes_deg_tmp){
  print(gene_tmp)
  for(cellClu_tmp in names(list_DEG_perClu)){
    df_cellType_clu_tmp = df[cellClu_annotation_kp == cellClu_tmp,]
    df_cellType_clu_tmp_back = df[cellClu_annotation_kp != cellClu_tmp,]
  
    df_cellType_clu_tmp$expr = df_cellType_clu_tmp[, gene_tmp]
    df_cellType_clu_tmp_back$expr = df_cellType_clu_tmp_back[, gene_tmp]
  
    df_cellType_clu_tmp$cell_cluster = cellClu_tmp
    df_cellType_clu_tmp_back$cell_cluster = cellClu_tmp
    
    ############# 
    
    pseudobulk_clu_tmp =  unlist(lapply(SampleIDs_kp,function(xx){
      mean(df_cellType_clu_tmp[which(df_cellType_clu_tmp$sample_id == xx),"expr"])
    }))
    
   
    df_pseudobulk_tmp = data.frame(expr = pseudobulk_clu_tmp,
                                   sample_id = SampleIDs_kp,
                                   sex = sex_perSampleID[SampleIDs_kp],
                                   cell_cluster = paste0("cell cluster ", cellClu_tmp))
    
    df_pseudobulk_tmp$sample_id = factor(df_pseudobulk_tmp$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
    df_pseudobulk_tmp$sex = factor(df_pseudobulk_tmp$sex, level = c("Male", "Female"))

    
    ############# 

    if(cellClu_tmp == names(list_DEG_perClu)[1]){
      df_  = df_pseudobulk_tmp
    }else{
      df_  = rbind(df_, df_pseudobulk_tmp)
      
    }
  }
  
  pdf(paste0("representative_plots/",gene_tmp,".pdf"),width=40,height=5)
  p_pseudobulk <- ggplot(df_, aes(x = sex, y = expr, fill =  sex)) + 
    theme_bw() +
    geom_violin(width=1,position=position_dodge(1)) +
    geom_boxplot(width=0.1, position=position_dodge(1)) +
    ggtitle(gene_tmp) +
    # title(gene_tmp)+
    theme(plot.title = element_text(face = "italic"),
          panel.grid = element_blank()
    ) + facet_wrap(~ cell_cluster, ncol = length(list_DEG_perClu))+geom_point()+
    scale_fill_manual(values=c( "#56B4E9","#E69F00" ))
  
  print(p_pseudobulk)
  
  dev.off()
  
  pdf(paste0("representative_plots/",gene_tmp,"_v2.pdf"),width=40,height=5)
  # print(p_pseudobulk + geom_point())
  p_pseudobulk <- ggplot(df_, aes(x = sex, y = expr, fill =  sex)) + 
    theme_bw() +
    # geom_violin(width=1,position=position_dodge(1)) +
    geom_boxplot(width=0.1, position=position_dodge(1)) +
    ggtitle(gene_tmp) +
    theme(plot.title = element_text(face = "italic"),
          panel.grid = element_blank()
    ) + facet_wrap(~ cell_cluster, ncol = length(list_DEG_perClu))+ geom_point()+
  scale_fill_manual(values=c( "#56B4E9","#E69F00" ))
  
  print(p_pseudobulk)
  dev.off()
  

}


## regular plot, for each cell cluster

for(cellClu_tmp in names(list_DEG_perClu)){
  print(cellClu_tmp)
  
  df_cellType_clu_tmp = df[cellClu_annotation_kp == cellClu_tmp,]
  df_cellType_clu_tmp_back = df[cellClu_annotation_kp != cellClu_tmp,]
  
  for(gene_tmp in genes_deg_tmp){
    dir.create(paste0("representative_plots/",gene_tmp))

    print(gene_tmp)
    
    plist = list()

    df_cellType_clu_tmp$expr = df_cellType_clu_tmp[, gene_tmp]
    df_cellType_clu_tmp_back$expr = df_cellType_clu_tmp_back[, gene_tmp]
    
    range_exp_geneTmp = range(df[,gene_tmp])
    
    pdf(paste0("clu", cellClu_tmp,"_",gene_tmp,".pdf"),width=40,height=40)
      
    
      p1 <- ggplot(df_cellType_clu_tmp, aes(x = x_location, y = y_location, color = expr)) +
        geom_point(size = .0005) +
        coord_fixed() +
        scale_color_gradient(low = "gray90", high = "blue",
                             trans = "sqrt", breaks = range_exp_geneTmp,
                             name = "expression level") +
        ggtitle(paste0(gene_tmp,", cluster ", cellClu_tmp)) + 
        theme_bw() +
        theme(plot.title = element_text(face = "italic"),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())+
        facet_wrap(~sex + sample_id, ncol = 7)+ theme_dark()
      
      p1_box <- ggplot(df_cellType_clu_tmp, aes(x = sex, y = expr, fill = sample_id)) + 
        ggtitle(paste0(gene_tmp,", cluster ", cellClu_tmp)) + 
        theme_bw() + 
        geom_violin(width=1,position=position_dodge(1)) +
        geom_boxplot(width=0.1, position=position_dodge(1)) +
        # geom_jitter()+
        # shape=16, position=position_jitter(0.2)
        theme(plot.title = element_text(face = "italic"),
              panel.grid = element_blank()
        )+ ylim(range_exp_geneTmp)+
        scale_fill_manual(values=c( "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9","#56B4E9","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00" ))
      
      
      
        p2 <- ggplot(df_cellType_clu_tmp_back, aes(x = x_location, y = y_location, color = expr)) +
        geom_point(size = .0005) +
        coord_fixed() +
        scale_color_gradient(low = "gray90", high = "blue",
                             trans = "sqrt", breaks = range_exp_geneTmp,
                             name = "expression level") +
        ggtitle(paste0(gene_tmp, ", outside cluster ", cellClu_tmp)) + 
        theme_bw() +
        theme(plot.title = element_text(face = "italic"),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())+
              facet_wrap(~sex + sample_id, ncol = 7)+ theme_dark()
            
      
      p2_box <- ggplot(df_cellType_clu_tmp_back, aes(x = sex, y = expr, fill = sample_id)) + 
        ggtitle(paste0(gene_tmp, ", outside cluster ", cellClu_tmp)) + 
        theme_bw() + 
        geom_violin(width=1,position=position_dodge(1)) +
        geom_boxplot(width=0.1, position=position_dodge(1)) +
        # geom_jitter()+
        # shape=16, position=position_jitter(0.2)
        theme(plot.title = element_text(face = "italic"),
              panel.grid = element_blank()
        ) + ylim(range_exp_geneTmp)+
        scale_fill_manual(values=c( "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9","#56B4E9","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00" ))
      
      
      ############# 
      
      pseudobulk_clu_tmp =  unlist(lapply(SampleIDs_kp,function(xx){
        mean(df_cellType_clu_tmp[which(df_cellType_clu_tmp$sample_id == xx),"expr"])
      }))
      
      pseudobulk_clu_tmp_back = unlist(lapply(SampleIDs_kp,function(xx){
        mean(df_cellType_clu_tmp_back[which(df_cellType_clu_tmp_back$sample_id == xx),"expr"])
                                }))
        
      df_pseudobulk_tmp = data.frame(expr = c(pseudobulk_clu_tmp, pseudobulk_clu_tmp_back),
                                     sample_id = rep(SampleIDs_kp, 2),
                                     sex = sex_perSampleID[rep(SampleIDs_kp, 2)],
                                     group = rep(c("within cluster","outside cluster"), each = length(SampleIDs_kp)))
      
      df_pseudobulk_tmp$sample_id = factor(df_pseudobulk_tmp$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
      df_pseudobulk_tmp$sex = factor(df_pseudobulk_tmp$sex, level = c("Male", "Female"))
      # table(df$sex)
      df_pseudobulk_tmp$group = factor(df_pseudobulk_tmp$group, level = c("within cluster","outside cluster"))
      
      p12_pseudobulk <- ggplot(df_pseudobulk_tmp, aes(x = sex, y = expr, fill = sex)) + 
        theme_bw() +
        geom_violin(width=1,position=position_dodge(1)) +
        geom_boxplot(width=0.1, position=position_dodge(1)) +
        geom_point() +
        theme(plot.title = element_text(face = "italic"),
              panel.grid = element_blank()
        ) + facet_wrap(~group, ncol = 2) +
        scale_fill_manual(values=c( "#56B4E9","#E69F00" ))
      
      ############# 
      print(ggarrange(p1, p1_box, p2, p2_box, p12_pseudobulk, heights = c(3, 0.7, 3, 0.7, 3),
                      ncol = 1, nrow = 5))
      
      dev.off()
      
    }

}


############# plot pseudobulk-level violin, for presentation
genes_deg_tmp = c("LYPD6B", "CDH13","CDH4")
gene_tmp = genes_deg_tmp[1]

cellClu_tmp = "33"

df_cellType_clu_tmp = df[cellClu_annotation_kp == cellClu_tmp,]
df_cellType_clu_tmp_back = df[cellClu_annotation_kp != cellClu_tmp,]

for(gene_tmp in genes_deg_tmp){
  dir.create(paste0("representative_plots/",gene_tmp))
  
  print(gene_tmp)
  
  plist = list()
  
  df_cellType_clu_tmp$expr = df_cellType_clu_tmp[, gene_tmp]
  df_cellType_clu_tmp_back$expr = df_cellType_clu_tmp_back[, gene_tmp]
  
  range_exp_geneTmp = range(df[,gene_tmp])
  
  ############# 
  
  pseudobulk_clu_tmp =  unlist(lapply(SampleIDs_kp,function(xx){
    mean(df_cellType_clu_tmp[which(df_cellType_clu_tmp$sample_id == xx),"expr"])
  }))
  
  pseudobulk_clu_tmp_back = unlist(lapply(SampleIDs_kp,function(xx){
    mean(df_cellType_clu_tmp_back[which(df_cellType_clu_tmp_back$sample_id == xx),"expr"])
  }))
  
  df_pseudobulk_tmp = data.frame(expr = c(pseudobulk_clu_tmp, pseudobulk_clu_tmp_back),
                                 sample_id = rep(SampleIDs_kp, 2),
                                 sex = sex_perSampleID[rep(SampleIDs_kp, 2)],
                                 group = rep(c("within cluster","outside of cluster"), each = length(SampleIDs_kp)))
  
  df_pseudobulk_tmp$sample_id = factor(df_pseudobulk_tmp$sample_id, level = SampleIDs_kp[order(sex_perSampleID[SampleIDs_kp])][nsample_kp:1])
  df_pseudobulk_tmp$sex = factor(df_pseudobulk_tmp$sex, level = c("Male", "Female"))
  # table(df$sex)
  df_pseudobulk_tmp$group = factor(df_pseudobulk_tmp$group, level = c("within cluster","outside of cluster"))
  
  p12_pseudobulk <- ggplot(df_pseudobulk_tmp, aes(x = sex, y = expr, fill = sex)) + 
    theme_bw() +
    geom_violin(width=1,position=position_dodge(1)) +
    geom_boxplot(width=0.1, position=position_dodge(1)) +
    geom_point() + ggtitle(gene_tmp)+
    theme(plot.title = element_text(face = "italic"),
          panel.grid = element_blank()
    ) + facet_wrap(~group, ncol = 2) +
    scale_fill_manual(values=c( "#00AFBB","#bb0c00" ))
  
  ############# 
  # print(ggarrange(p1, p1_box, p2, p2_box, p12_pseudobulk, heights = c(3, 0.7, 3, 0.7, 3),
                  # ncol = 1, nrow = 5))
  pdf(paste0("clu", cellClu_tmp,"_",gene_tmp,".pdf"),width=5,height=5)
  
  print(p12_pseudobulk)
  dev.off()
  
}












