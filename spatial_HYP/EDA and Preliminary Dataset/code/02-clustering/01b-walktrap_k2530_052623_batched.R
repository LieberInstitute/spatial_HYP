library(data.table)
library(Biostrings)
library(SpatialExperiment)
library(ggspavis)
library(scater) # addPerCellQC
library(scran)
library(igraph)

hyp2 <- readRDS("hyp_umi600_gene450_chrm35_lognorm.RDS")
hyp2.vgs.list <- readRDS("nnSVG_HVG_and_combo_featurelists.RDS")
glists2 <- readRDS("SNNgraphs_svgsAndHvgs_k25_k30_WITHOUTcovariates_052623.RDS")
pcalist <- readRDS("nnsvg_HVG_and_combos_PCA_NOcovariates.RDS")

taskid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k <- ceiling(taskid/7)

if (taskid>7){taskid<-taskid-7}
### This block still wouldn't run locally even with better mem mgmt. Save relevant files and run through JHPCE. Each feature set is trying to use about 20 GB for the k=25 (so ~140 GB total for k25), and probably 25% more for the k=30 (175 GB). glists2[[]] is just the k=25 and k=30 SNNgraphs so run the commented out block below accordingly. 


tmpwalk <- igraph::cluster_walktrap(glists2[[k]][[taskid]])
walkclusts.tmpout <- tmpwalk$membership
     
 
### split run only, continued
 
saveRDS(walkclusts.tmpout,paste0("out/",names(glists2)[k],"_",names(pcalist)[taskid],"_walktrap_clusters_nocovar.RDS"))

