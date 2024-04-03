library(data.table)

tired <- fread("~/Desktop/tiredduringdep.tsv")

wgain <- fread("~/Desktop/weightgainworst.tsv")

wloss <- fread("~/Desktop/weightlossduringworstdep.tsv")

slp <- fread("~/Desktop/depdeltasleep.tsv")

gw <- list(tired,wgain,wloss,slp)
names(gw) <- c("tiredness","wgain","wloss","slpchange")
rm(tired,wgain,wloss,slp)

gc(full=T)


gw <- lapply(gw,FUN=function(x){
	x[,one_allele:=gsub(variant,pattern="^.*:.*:(.*):.*$",replacement="\\1")]
	x[,two_allele:=gsub(variant,pattern="^.*:.*:.*:(.*)$",replacement="\\1")]
	x[,majal:=""]
	x[one_allele==minor_allele,majal:=two_allele]
	x[two_allele==minor_allele,majal:=one_allele]})


gw <- lapply(gw,FUN=function(x){x[,one_allele:=NULL]
	x[,two_allele:=NULL]})

gc(full=T)


gw <- lapply(gw,FUN=function(x){x[,variant2:=gsub(variant,pattern="^(.*):(.*):(.*):(.*)",replacement=paste0("\\1",":","\\2","_","\\3","_","\\4"))]})

lut <- fread("~/Desktop/full_variant_qc_metrics.txt",select=c(5,6))
setnames(lut,c("rsid","crap"))

# max file size is 75MB for webtwas; get there by cutting to unadj p<0.25
gw <- lapply(gw,FUN=function(x){
	y <- merge.data.table(x,lut,by.x="variant2",by.y="crap")
	y <- y[pval<0.25]
	y <- y[,.(rsid,minor_allele,majal,beta,pval)]
	return(y)})

rm(lut)
gc(full=T)

mapply(X=names(gw),Y=gw,FUN=function(X,Y){fwrite(Y,file=paste0("~/Desktop/",X,"mod.tsv.gz"),row.names=F,col.names=T,sep='\t',quote=F)})