setwd("/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/")
library(data.table)
library(zellkonverter)
library(SingleCellExperiment)

cellmeta <- fread("raw-data/ABA_Whole_Mouse_Brain_Yao23_HYP/cell_metadata_with_cluster_annotation.csv")
hypv2 <- zellkonverter::readH5AD("raw-data/ABA_Whole_Mouse_Brain_Yao23_HYP/WMB-10Xv2-HY-raw.h5ad")
hypv3 <- zellkonverter::readH5AD("raw-data/ABA_Whole_Mouse_Brain_Yao23_HYP/WMB-10Xv3-HY-raw.h5ad")

tmpcdv2 <- as.data.table(colData(hypv2),keep.rownames=T)
cellmeta2 <- cellmeta[cell_label %in% colnames(hypv2)]
tmpcdv2 <- merge.data.table(tmpcdv2,cellmeta2,by.x="rn",by.y="cell_label")
tmpcdv2 <- DataFrame(tmpcdv2)
rownames(tmpcdv2) <- tmpcdv2$rn
tmpcdv2 <- tmpcdv2[colnames(hypv2),]
colData(hypv2) <- tmpcdv2

tmpcdv3 <- as.data.table(colData(hypv3),keep.rownames=T)
cellmeta3 <- cellmeta[cell_label %in% colnames(hypv3)]
tmpcdv3 <- merge.data.table(tmpcdv3,cellmeta3,by.x="rn",by.y="cell_label")
tmpcdv3 <- DataFrame(tmpcdv3)
rownames(tmpcdv3) <- tmpcdv3$rn
tmpcdv3 <- tmpcdv3[colnames(hypv3),]
colData(hypv3) <- tmpcdv3

## since these are all uniformly processed for allen atlas, should be able to join them up...
bothhyp <- cbind(hypv2,hypv3)
assayNames(bothhyp) <- "counts"
rm(tmpcdv2,tmpcdv3,cellmeta2,cellmeta3,cellmeta,hypv2,hypv3)
gc(full=T)

colData(bothhyp) <- colData(bothhyp)[,c(1:3,8,10,12,15,16,20,21)]
colnames(colData(bothhyp))
bothhyp$donor_sex <- as.factor(bothhyp$donor_sex)
rownames(bothhyp) <- rowData(bothhyp)$gene_symbol

# drop genes with two or more rows
drops <- as.data.table(rownames(bothhyp))[,.N,by="V1"]
keeprows <- drops[N==1,V1]

bothhyp <- bothhyp[keeprows,]
stopifnot(identical(rownames(bothhyp),rowData(bothhyp)$gene_symbol)&nrow(bothhyp)==length(unique(rownames(bothhyp))))

rm(keeprows,drops)
gc(full=T)

### fix supercluster/subclass labels to match those from registration
# Fix variable names for compatibility
bothhyp$subclass <- gsub(bothhyp$subclass,pattern=" ",replacement="_")
bothhyp$subclass <- gsub(bothhyp$subclass,pattern="-",replacement="_")
bothhyp$subclass <- gsub(bothhyp$subclass,pattern="/",replacement="_")
bothhyp$subclass <- paste0("x",bothhyp$subclass)
bothhyp$subclass <- as.factor(bothhyp$subclass)

bothhyp$supertype <- gsub(bothhyp$supertype,pattern=" ",replacement="_")
bothhyp$supertype <- gsub(bothhyp$supertype,pattern="-",replacement="_")
bothhyp$supertype <- gsub(bothhyp$supertype,pattern="/",replacement="_")
bothhyp$supertype <- paste0("x",bothhyp$supertype)
bothhyp$supertype <- as.factor(bothhyp$supertype)

## save object (now a Singlecellexperiment by way of zellkonverter) for processing on local machine, where dreamlet works but zellkonverter doesn't.

saveRDS(bothhyp,file="raw-data/ABA_Whole_Mouse_Brain_Yao23_HYP/bothdatasets_raw_asSCE.RDS")

sessionInfo()
