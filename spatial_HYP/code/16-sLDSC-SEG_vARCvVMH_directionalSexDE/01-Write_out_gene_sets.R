
setwd("/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/spatial_HYP")
library(data.table)

## make a table with the cluster annots
anno <- as.data.table(as.data.frame(matrix(nrow=15,ncol=1)))
setnames(anno,c("clusid"))
anno[,clusid:=paste0("X",c(1:15))]
anno[clusid=="X7",anno:="VMH.1"]
anno[clusid=="X12",anno:="VMH.2"]
anno[clusid=="X4",anno:="ARC.1"]
anno[clusid=="X6",anno:="ARC.2"]
anno[clusid=="X3",anno:="OT.1"]
anno[clusid=="X10",anno:="OT.2"]
anno[clusid=="X5",anno:="OT.3"]
anno[clusid=="X1",anno:="GABA.1"]
anno[clusid=="X2",anno:="PeriVN"]
anno[clusid=="X9",anno:="Vascular"]
anno[clusid=="X8",anno:="SON"]
anno[clusid=="X13",anno:="PortalVasc"]
anno[clusid=="X14",anno:="Astro"]
anno[clusid=="X15",anno:="GABA.2"]
anno[clusid=="X11",anno:="DROP"] # sample specific cluster

sexde.sing <- fread("processed-data/09-Sex DE/01b-voomLmFit_nnsvg10-HmnylmbNA-BS-15-singleClusts_ARC2fixed.txt")
sexde.jt <- fread("processed-data/09-Sex DE/01c-voomLmFit_nnsvg10-HmnylmbNA-BS-15-VMHARCclpsd.txt")

sexde.jt <- sexde.jt[assay %in% c("ARC","VMH")]
sexde.jt <- sexde.jt[,assay:=paste0("v",assay)]
sexde <- as.data.table(rbind(sexde.jt,sexde.sing))

setnames(sexde,c("gene_id","gene_name"),c("ensembl","gene"))
allg <- as.data.table(unique(sexde[,.(ensembl,gene)]))

sexde <- merge.data.table(sexde,anno,by.x="assay",by.y="clusid",all.x=T)
sexde[assay %in% c("vARC","vVMH"),anno:=assay]
sexde <- sexde[assay %in% c("vARC","vVMH","VMH.1","VMH.2","ARC.1","ARC.2")]

setnames(sexde,"t","t_stat")

rm(sexde.jt,sexde.sing)

allg <- as.data.table(unique(sexde[,.(ensembl,gene)]))

for (ann in unique(sexde$anno)){
    tmpdat <- copy(sexde[anno==ann])
    setorderv(tmpdat,"t_stat",-1) # descending, ie male up
    tmpdat <- tmpdat[c(1:ceiling(0.05*nrow(tmpdat)))]
    

    allg[,newcol:=0]
    allg[ensembl %in% tmpdat$ensembl,newcol:=1]
    setnames(allg,"newcol",paste0(gsub(ann,pattern="\\.",replacement=""),"_maleup"))
    rm(tmpdat)

    tmpdat <- copy(sexde[anno==ann])
    setorderv(tmpdat,"t_stat",1) # ascending, ie fem up
    tmpdat <- tmpdat[c(1:ceiling(0.05*nrow(tmpdat)))]

    allg[,newcol:=0]
    allg[ensembl %in% tmpdat$ensembl,newcol:=1]
    setnames(allg,"newcol",paste0(gsub(ann,pattern="\\.",replacement=""),"_femup"))
    rm(tmpdat)
}

hg19g <- fread("/dcs04/lieber/shared/statsgen/LDSC/base/gene_meta_hg19.txt")
setnames(hg19g,c("ensg","ensgver","chr","start","end","gene_name","pctgc","gtype","entrez"))

for (n in names(allg)[3:ncol(allg)]){
    tmp <- copy(hg19g[ensg %in% allg[get(n)==1,ensembl]])[,.(ensg,chr,start,end)]
    tmp[,start:=ifelse(start-100000<0,0,start-100000)]
    tmp[,end:=end+100000]
    tmp[,newchr:=paste0("chr",chr)]
    tmp[,chr:=NULL]
    setnames(tmp,"newchr","chr")
    tmp <- tmp[!(chr %in% paste0("chr",c("Y","X","MT")))]
    stopifnot(nrow(tmp[start>end])==0)
    tmp <- tmp[,.(chr,start,end)]
    fwrite(tmp,paste0("code/16-sLDSC-SEG_test/bedfiles/",n,".bed"),sep='\t',quote=F,row.names=F,col.names=F)
}
