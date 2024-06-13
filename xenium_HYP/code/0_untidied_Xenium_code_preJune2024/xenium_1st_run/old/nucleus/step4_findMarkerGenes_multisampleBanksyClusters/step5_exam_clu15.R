# ml conda_R/4.3
# R

###########
lambda <- 0
lambdaName="0"

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

################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)

se_seurat  <- as.Seurat(se) 
Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
# se_seurat = FindAllMarkers(se_seurat)

# saveRDS(se_seurat, file = paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, "_nucleus.rds"))

# se_seurat = readRDS(paste0(dir_data_processed_banksy, "se_seurat_markerGenes_multiSampleBanksy_lambda", lambdaName, "_nucleus.rds"))
######################################################################### 
### plot the results - heatmap for marker genes
######################################################################### 
library(dplyr)
genes_interest = c("gad1", "gad2", "slc32a1", "slc17a6",  "slc17a7")
genes_interest = toupper(genes_interest)

######## process data into seurat format
counts_clu15 = assays(se)$counts[,colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]=="15"]

summary(counts_clu15["GAD1",])
summary(counts_clu15["SLC17A6",])

######## 
obj = CreateSeuratObject(assays(se)$counts)
obj = NormalizeData(obj)
obj = ScaleData(obj)
######## subset of cluster-15 cells
obj_clu15 = obj[,Idents(se_seurat)=="15"]
data_clu15 =  as.matrix(obj_clu15@assays$RNA@data)
dim(data_clu15)
# [1]  366 4132
a=assays(se)$counts["GAD1",Idents(se_seurat)=="15"]==0
summary(data_clu15["GAD1",a])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       0       0       0       0 
summary(data_clu15["GAD1",])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   4.687   3.330   5.649   7.355 


counts_Gad_clu15 = assays(se)$counts["GAD1",Idents(se_seurat)=="15"]
summary(counts_Gad_clu15)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   1.000   4.075   6.000  55.000

######## check the overall expressing level of the genes in clu 15
summary(data_clu15["GAD1",])
summary(data_clu15["GAD2",])
summary(data_clu15["SLC17A6",])
summary(data_clu15["SLC17A7",])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   4.687   3.330   5.649   7.355 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.3297  0.0000  6.2166 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   3.606   4.841   3.908   5.398   6.909 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.00000 0.02209 0.00000 5.77945

# GAD1，SLC17A6 both express highly expressing in clu15

######## check if there's a cell that's high in both sla and gad
# cells_highGAD1 = colnames(data_clu15)[data_clu15["GAD1",]>6]
dim(data_clu15)
# [1]  366 4132
# 
# cells_highSLC17A6 = colnames(data_clu15)[data_clu15["SLC17A6",]>6]
sum(df_[,1]==0&df_[,2]==0)
sum(df_[,1]>0&df_[,2]>0)
sum(df_[,1]>0&df_[,2]==0)
sum(df_[,1]==0&df_[,2]>0)
# 412
# 1] 2070
# [1] 516
# [1] 1134

######### scatter plot GAD vs SL

df_=data.frame(GAD1 =data_clu15["GAD1",],
               SLC17A6 =data_clu15["SLC17A6",] )
pdf("scatter_GAD1_SLC17A6_cl15.pdf")
# plot(data_clu15["GAD1",], data_clu15["SLC17A6",],
     # xlab = "GAD1 expression level", ylab = "SLC17A6 expression level")
ggplot(df_, aes(x=GAD1, y=SLC17A6)) + 
  geom_point(size=.01,aes(alpha=.001))+xlab("GAD1 expression level")+ylab("GAD1 expression level")
dev.off()

######## SP plot, plot those cells in the space
df_=data.frame(GAD1 =data_clu15["GAD1",],
               SLC17A6 =data_clu15["SLC17A6",] )

library(ggpubr)

clu_tmp = "15"

list_cells_which = list(bothHigh = df_[,1]>0&df_[,2]>0,
                   SLC17A6highOnly= df_[,1]==0&df_[,2]>0,
                   GAD1highOnly = df_[,1]>0&df_[,2]==0 )

dir_out_clu_tmp = paste0("cells_exam_clu", clu_tmp)
dir.create(dir_out_clu_tmp)
range_ = range(as.matrix(df_))

for(name_grpCell in names(list_cells_which)){
  
  spe=se[,colnames(data_clu15)[list_cells_which[[name_grpCell]]]]
  
  list_p = list()
  png(paste0(dir_out_clu_tmp, "/", name_grpCell, ".png"), res=200,
      width=30, height=15, units="in")
  
  for(gene_tmp in c("GAD1", "SLC17A6")){
    
    df <- as.data.frame(
      cbind(spatialCoords(spe), 
            expr = logcounts(spe)[gene_tmp, ]))
    df$sample_id=colData(spe)$sample_id
    
    df$sex=sex_perSampleID[as.character(df$sample_id)]
    df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
    
    xlim_ = range(df$x_location)
    ylim_ = range(df$y_location)
    
    list_p[[gene_tmp]] <- ggplot(df, 
                aes(x = x_location, y = y_location, color = expr)) +
      geom_point(size = .0005) +
      coord_fixed()+
      ggtitle(paste0(gene_tmp,", cells ",name_grpCell) )+
      # theme_bw() +
      xlim(xlim_)+ylim(ylim_)+
      theme(plot.title = element_text(size =55), 
            panel.background = element_rect(fill = "gray24",
                                            colour = "gray24"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())+
      facet_wrap(~sex + sample_id, ncol = 3) +
      scale_color_gradient(low = "white", high = "blue",
                           trans = "sqrt",
                           breaks = range_, 
                           limits = range_,
                           name = "expression level") 
    
    # print(p)
    
  }
  
  print(ggarrange(list_p[[1]], list_p[[2]], heights = c(2, 2),
                  ncol = 1, nrow = 2))
  dev.off()
  
}
# 