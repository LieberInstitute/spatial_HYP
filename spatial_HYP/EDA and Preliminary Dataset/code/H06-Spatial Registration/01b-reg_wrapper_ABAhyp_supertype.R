
#### sets tab autocompletion for directories to begin from the directory containing the .Rproj session, instead of the script's directory if different

library(data.table)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(SpatialExperiment)
library(spatialLIBD)
library(zellkonverter)
library(magrittr)
library(BiocParallel)

setDTthreads(7)
### data sources
# https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/WMB-10X/20230630/views/
# https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-HY-log2.h5ad
# https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/WMB-10X/20230630/views/
# https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv2/20230630/WMB-10Xv2-HY-log2.h5ad
# https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-HY-raw.h5ad
# https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv2/20230630/WMB-10Xv2-HY-raw.h5ad


### load counts

cellmeta <- fread("ABA_Whole_Mouse_Brain_Yao23_HYP/cell_metadata_with_cluster_annotation.csv")
hypv2 <- zellkonverter::readH5AD("ABA_Whole_Mouse_Brain_Yao23_HYP/WMB-10Xv2-HY-raw.h5ad")
hypv3 <- zellkonverter::readH5AD("ABA_Whole_Mouse_Brain_Yao23_HYP/WMB-10Xv3-HY-raw.h5ad")

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

# hypv2log <- zellkonverter::readH5AD("ABA_Whole_Mouse_Brain_Yao23_HYP/WMB-10Xv2-HY-log2.h5ad")
# hypv3log <- zellkonverter::readH5AD("ABA_Whole_Mouse_Brain_Yao23_HYP/WMB-10Xv3-HY-log2.h5ad")
#
# stopifnot(identical(c(colnames(hypv2log),colnames(hypv3log)),colnames(bothhyp)))
# hyplogs <- cbind(hypv2log,hypv3log)
#
# assays(bothhyp)$logcounts <- assays(hyplogs)[[1]]
# rm(hypv2log,hypv3log,hyplogs)
# gc(full=T)


# Filter to abundant enough cells then prep for spatial reg

tab <- as.data.table(table(bothhyp$library_label.x,bothhyp$supertype))
tab <- dcast(tab,V2~V1,value.var = "N")

# subset to cell subclasses or supertypes with at least 5 cells in at least 5 samples
tab[,sum:=rowSums(tab[,.SD>=5,.SDcols=c(2:ncol(tab))])]
# subsetting for # of sample•(supertypes or subclasses)
tab <- tab[sum>4]

clusts <- as.list(tab$V2)

# get column names of single cells from samples contributing at least 5 cells of the current (subclasses or supertypes)
keepcols <- unique(unlist(sapply(clusts,FUN=function(x){
    keepsamps <- names(tab)[which(tab[V2==x,.SD>5,.SDcols=c(2:ncol(tab)-1)])]
    colnames(bothhyp[,bothhyp$supertype==x&bothhyp$library_label.x %in% keepsamps])
},simplify = T)))
rm(clusts)

bothhyp <- bothhyp[,keepcols]
gc(full=T)

# Fix variable names for compatibility
bothhyp$supertype <- gsub(bothhyp$supertype,pattern=" ",replacement="_")
bothhyp$supertype <- gsub(bothhyp$supertype,pattern="-",replacement="_")
bothhyp$supertype <- gsub(bothhyp$supertype,pattern="/",replacement="_")
bothhyp$supertype <- paste0("x",bothhyp$supertype)
bothhyp$supertype <- as.factor(bothhyp$supertype)

# fix genotypes to be simple (there's SNAP25 cre reporter and Ai14-tdt wt) so wew can also use them as a covariate
bothhyp$donor_genotype[which(bothhyp$donor_genotype=="Snap25-IRES2-Cre/wt;Ai14(RCL-tdT)/wt")] <- "snapTdt"
bothhyp$donor_genotype[which(bothhyp$donor_genotype=="Ai14(RCL-tdT)/wt")] <- "wtTdt"
## wait nvm, we don't want genotype (which was used to label cells for sorting and collection across several different genotypes across the larger dataset) as a covariate--its going to be confounded with cell type in this case (all cells or neurons only)
# bothhyp$donor_genotype <- as.factor(bothhyp$donor_genotype)

# set library method (v2 or v3) as factor (which needs a leading character not digit) so we can control for it
bothhyp$library_method[which(bothhyp$library_method=="10Xv2")] <- "v2"
bothhyp$library_method[which(bothhyp$library_method=="10Xv3")] <- "v3"
bothhyp$library_method <- as.factor(bothhyp$library_method)

# make rownamnes capitalized symbols so we can cross-map to human
rowData(bothhyp)$gene_symbol <- toupper(rowData(bothhyp)$gene_symbol)
rownames(bothhyp) <- rowData(bothhyp)$gene_symbol
# drop genes with two or more rows
drops <- as.data.table(rownames(bothhyp))[,.N,by="V1"]
keeprows <- drops[N==1,V1]

bothhyp <- bothhyp[keeprows,]
stopifnot(identical(rownames(bothhyp),rowData(bothhyp)$gene_symbol)&nrow(bothhyp)==length(unique(rownames(bothhyp))))

message("Input dataset ready. Starting spatialLIBD::registration_wrapper...")

abahypres <- spatialLIBD::registration_wrapper(sce = bothhyp,var_registration = "supertype",var_sample_id = "library_label.x",gene_ensembl = "gene_symbol",gene_name = "gene_symbol",min_ncells = 5)

saveRDS(abahypres,"out/ABAhyp_mincells5-n5_supertype_modelingres.RDS")


