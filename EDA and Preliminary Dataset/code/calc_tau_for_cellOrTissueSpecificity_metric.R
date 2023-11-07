### original paper: Yanai et al 2005 Bioinformatics
### follow-up showing this is the most robust/least biased metric for tissue specificity, still, in 2017:  https://academic.taup.com/bib/article/18/2/205/2562739#119555146
### use of this statistic on (pseudobulked) cell types from scRNAseq: Allen Brain Adult Mouse scRNAseq Yao van Velthoven..Zeng 23 preprint

calc_tau <- function(table, genes.in.rows=T,id.var.for.rejoining.values="rn"){
    if (genes.in.rows==FALSE){
        table <- as.data.table(t(table),keep.rownames=T)
    }
    if (identical(rownames(table),seq(1,nrow(table),1))){
        working <- as.data.table(table)
    }
    else {working <- as.data.table(table,keep.rownames=T)}
    usecols <- names(working)[!(names(working) %in% c("rn",id.var.for.rejoining.values))]
    working[,maxv:=apply(.SD,1,max),.SDcols=usecols]
    set(working,j=usecols,value=working[,..usecols]/working[,maxv])
    working[,tau:=apply(.SD,1,FUN=function(x){sum(1-x)/(length(x)-1)}),.SDcols=usecols]
    returnnames <- names(working)[names(working) %in% c("rn",id.var.for.rejoining.values,"tau")]
    return(working[,..returnnames])
}

