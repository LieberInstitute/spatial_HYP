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
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))
se = spe_joint

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
              "br1735a","br1735b","br6588")
sex_perSampleID = c("M", "F", "F", "M", "M", "F")
names(sex_perSampleID) = sampleID_unique
# sex_perSampleID[sampleID_unique]
########################################################################
################################# find marker genes
######################################################################### 
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)
cellClu = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

se_seurat  <- as.Seurat(se) 
se_seurat@meta.data$sex = sex_perSampleID[colData(se)$sample_id]
sex = sex_perSampleID[colData(se)$sample_id]

####### VHM DEG
spe = se_seurat[,cellClu=="4"]


# Idents(se_vhm) <- se_vhm@meta.data$sex
# se_vhm = FindAllMarkers(se_vhm)
# saveRDS(se_vhm, file = paste0(dir_data_processed_banksy, 
#                               "se_vhm_seurat_sexDEG_multiSampleBanksy_lambda", 
#                               lambdaName, 
#                               ".rds"))
# se_vhm["ADARB1",]

####### ARC DEG
spe = se[,cellClu=="3"]

pdf("try.pdf",width=30,height=30)
gene_tmp = "RXRG"
df <- as.data.frame(
    cbind(spatialCoords(spe), 
          expr = logcounts(spe)[gene_tmp, ]))
df$sample_id=colData(spe)$sample_id

df$sex=sex_perSampleID[as.character(df$sample_id)]
df$sex_sampleID=paste0(df$sex,"_",df$sample_id)
table(df$sex,df$sample_id) ##WTF???

ggplot(df, aes(x = x_location, y = y_location, color = expr)) + 
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
dev.off()


# Idents(se_arc) <- se_arc@meta.data$sex
# se_arc = FindAllMarkers(se_arc)

# saveRDS(se_arc, file = paste0(dir_data_processed_banksy, 
                              # "se_arc_seurat_sexDEG_multiSampleBanksy_lambda", 
                              # lambdaName, ".rds"))


# ADARB1
se_arc["ADARB1",]
######################################################################### 
### plot the results - heatmap for marker genes
######################################################################### 
library(dplyr)
######################################## 
###################### plot vhm
######################################## 
cellClu_used = "4"

pdf(paste0("DE_gender_vhm_multisampleBanksy_lambda_",lambdaName,".pdf"), width = 20, height = 20)
  ######## format se into seurat format
  obj = CreateSeuratObject(assays(se)$counts[, cellClu == cellClu_used], 
                           project = paste0("marker genes"))
  obj = NormalizeData(obj)
  obj = ScaleData(obj)
  Idents(obj) = sex[cellClu == cellClu_used]
    
  ######## 
  se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, 
                                     "se_vhm_seurat_sexDEG_multiSampleBanksy_lambda",
                                     lambdaName,".rds"))
  
  se_seurat_markers %>%
    group_by(cluster) %>%
    # dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  

  set.seed(111)#plot a subset of cells to reduce running time
  rd_cells = sample(1:ncol(obj), 10000)
  print(DoHeatmap(obj[top10$gene, rd_cells] , 
                  features = top10$gene) + NoLegend())

  # print(DoHeatmap(obj[, ] , 
                  # features = top10$gene) + NoLegend())
  
dev.off()

######################################## 
###################### plot arc
######################################## 
cellClu_used = "3"

pdf(paste0("DE_gender_arc_multisampleBanksy_lambda_",lambdaName,".pdf"), width = 20, height = 20)
######## format se into seurat format
obj = CreateSeuratObject(assays(se)$counts[, cellClu == cellClu_used], 
                         project = paste0("marker genes"))
obj = NormalizeData(obj)
obj = ScaleData(obj)
Idents(obj) = sex[ cellClu == cellClu_used]

######## 
se_seurat_markers = readRDS(paste0(dir_data_processed_banksy, 
                                   "se_arc_seurat_sexDEG_multiSampleBanksy_lambda",
                                   lambdaName,".rds"))

se_seurat_markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


set.seed(111)#plot a subset of cells to reduce running time
rd_cells = sample(1:ncol(obj), 10000)
print(DoHeatmap(obj[top10$gene, rd_cells] , features = top10$gene) + NoLegend())


dev.off()

