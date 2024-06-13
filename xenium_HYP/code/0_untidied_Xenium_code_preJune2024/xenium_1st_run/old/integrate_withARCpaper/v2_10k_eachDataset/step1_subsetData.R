dir_Data_arcPaper = "/dcs04/hansen/data/ywang/xenium/data_paperARC/"
########## load counts matrix
counts_arcPaper = Matrix::readMM(paste0(dir_Data_arcPaper,"GSE169109_matrix.mtx"))
# [1]  32738 135731
counts_arcPaper[1:2,1:2]
set.seed(111)
rd_=sample(1:ncol(counts_arcPaper),10000)
dim(counts_arcPaper)
# [1]  32738 135731

# saveRDS(counts_arcPaper, file = paste0(dir_Data_arcPaper,"counts_arcPaper.rds"))
counts_arcPaper_sub =counts_arcPaper[,rd_]
saveRDS(counts_arcPaper_sub, file = paste0(dir_Data_arcPaper,"counts_arcPaper_sub.rds"))
# saveRDS(counts_arcPaper_sub, file = paste0(dir_Data_arcPaper,"counts_arcPaper_sub.rds"))
# counts_arcPaper_sub = readRDS("counts_arcPaper_sub.rds")

########## load barcode
barcode_arcPaper  = read.delim(paste0(dir_Data_arcPaper,"GSE169109_barcodes.tsv"),header = F)
dim(barcode_arcPaper)
# [1] 135731      1

##########  load features
# GSE169109_features.tsv
feature_arcPaper  = read.delim(paste0(dir_Data_arcPaper,"GSE169109_features.tsv"),header = F)
dim(feature_arcPaper)
# [1] 32738     3
