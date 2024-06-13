dir_Data_arcPaper = "/dcs04/hansen/data/ywang/xenium/data_paperARC/"
setwd(dir_Data_arcPaper)
########## load counts matrix

set.seed(111)
rd_=sample(1:(135731),10000)

########## 
# saveRDS(counts_arcPaper_sub, file = paste0(dir_Data_arcPaper,"counts_arcPaper_sub.rds"))
counts_arcPaper_sub = readRDS("counts_arcPaper_sub.rds")

########## load barcode
barcode_arcPaper  = read.delim(paste0(dir_Data_arcPaper,"GSE169109_barcodes.tsv"),header = F)
dim(barcode_arcPaper)
# [1] 135731      1
barcode_arcPaper_sub = barcode_arcPaper[rd_,]
##########  load features
# GSE169109_features.tsv
feature_arcPaper  = read.delim(paste0(dir_Data_arcPaper,"GSE169109_features.tsv"),header = F)
dim(feature_arcPaper)
# [1] 32738     3

######### assign rownames and colnames to the count matrix
colnames(counts_arcPaper_sub) = barcode_arcPaper_sub
rownames(counts_arcPaper_sub) = toupper(feature_arcPaper[,2])
saveRDS(counts_arcPaper_sub, file = paste0(dir_Data_arcPaper,"counts_arcPaper_sub_formatted.rds"))

######### format into seurat obj
library(Seurat)
obj_arcPaper = CreateSeuratObject(counts_arcPaper_sub, project = "ARCpaper")
obj_arcPaper = NormalizeData(obj_arcPaper)
obj_arcPaper = FindVariableFeatures(obj_arcPaper)
saveRDS(obj_arcPaper, file=paste0(dir_Data_arcPaper,"obj_arcPaper_sub.rds"))



