######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 1

###########
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)

dir_data_processed = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
# dir_out = "/dcs04/hansen/data/ywang/xenium_newData/out/"
# setwd(dir_out)
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"

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

###########
# sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
###########
sample_names_old = readRDS(paste0("/dcs04/hansen/data/ywang/xenium/data_processed_banksy/", "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda1_13samples.rds"))

# spe_joint_2runs_banksy_lambda0.rds


##################################################################
##################################################################
########### find genes that are annotated
##################################################################
##################################################################

#######################################################
geneInfo_ = read.csv(paste0("/dcs04/hansen/data/ywang/xenium/data_processed/","data/Xenium_HYP_Candidate_Genes.csv"))
head(geneInfo_)

geneInfo_[,2]

genes_ARC = geneInfo_[grep("ARC",geneInfo_[,2]),]
genes_ARC_marker1 = genes_ARC[grep("ARC Enriched",genes_ARC[,2]),]
genes_ARC_marker2 = genes_ARC[grep("Known ARC marker",genes_ARC[,2]),]
genes_ARC_marker3 = genes_ARC[grep("ARC enriched",genes_ARC[,2]),]
genes_ARC_marker = rbind(genes_ARC_marker1, genes_ARC_marker2, genes_ARC_marker3)
genes_ARC_marker = intersect( genes_ARC_marker[,1], rownames(spe_joint))


genes_VMH = geneInfo_[grep("VMH",geneInfo_[,2]),]
genes_VMH_marker = genes_VMH[grep("DEG: VMH vs not",genes_VMH[,2]),]
genes_VMH_marker = intersect( genes_VMH_marker[,1], rownames(spe_joint))
# > genes_VMH_marker
# [1] "ANKRD34B"  "TLL2"      "CDH13"     "LINC01615" "CELF2"     "HRH3"     
# [7] "FRMPD2"    "GDF10"     "GLRA3"     "ICAM5"     "LBHD2"     "MYZAP"    
# [13] "BMP8B"     "DDN"       "NR2F2"     "RBFOX1"    "GABRA5"    "ELK1"     
# [19] "LAMP5"     "TMEM233"   "PCDH19"    "DOC2A"     "KRT8"      "KRT18"    
# > 
#   > genes_ARC_marker
# [1] "GLRA2"   "AGRP"    "SSTR1"   "GSX1"    "SLC10A4" "NPY2R"   "ENOX2"  
# [8] "NPY1R"   "CYP26B1" "MC3R"   


######################################################################### 
## find domains
######################################################################### 
library(mgcv)
x1 = colData(spe_joint)$x_position
x2 = colData(spe_joint)$y_position
## VMH
normcounts = assays(spe_joint)$normcounts

### ARC
exp_ARC = colSums(normcounts[genes_ARC_marker,])
length(exp_ARC)
# 911435
# gam: https://noamross.github.io/gams-in-r-course/chapter3

df_ = data.frame(x1=x1,x2=x2,exp_ARC=exp_ARC)
for(sample in sample_names[1:5]){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  fit_gam_arc = gam(exp_ARC ~ s(x1, x2),data = df_[which_sample,], # ~5min to fit
                    method = "REML")
  saveRDS(fit_gam_arc, file = paste0("fit_gam_arc_sample",sample,".rds") )
  pred_arc = predict.gam(fit_gam_arc) # ~1min to predict
  saveRDS(pred_arc, file = paste0("pred_gam_arc_sample",sample,".rds") )
  
}




pred_arc_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  pred_arc  =readRDS( paste0("pred_gam_arc_sample",sample,".rds") )
  pred_arc_c[which_sample] = pred_arc
}
set.seed(111)
domain_if_ARC = kmeans(pred_arc_c, 2)$cluster
colData(spe_joint)$domain_if_ARC = domain_if_ARC
saveRDS(domain_if_ARC, file = paste0("domain_if_ARC.rds") )



domain_if_ARC_v2= array(dim=ncol(spe_joint))
for(sample in sample_names){
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  pred_arc  =readRDS( paste0("pred_gam_arc_sample",sample,".rds") )
  set.seed(111)
  domain_if_ARC_sample_tmp =  kmeans(pred_arc, 2)$cluster
  domain_if_ARC_v2[which_sample] = domain_if_ARC_sample_tmp
}
# domain_if_VMH_v2 = kmeans(pred_vmh_c_v2, 2)$cluster
colData(spe_joint)$domain_if_ARC_v2 = domain_if_ARC_v2
saveRDS(domain_if_ARC_v2, file = paste0("domain_if_ARC_v2.rds") )


## VMH
exp_VMH = colSums(normcounts[genes_VMH_marker,])
# colData(spe)$domain_if_VMH
df_$exp_VMH = exp_VMH
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  fit_gam_vmh = gam(exp_VMH ~ s(x1, x2),data = df_[which_sample,], # ~5min to fit
                    method = "REML")
  saveRDS(fit_gam_vmh, file = paste0("fit_gam_vmh_sample",sample,".rds") )
  pred_vmh = predict.gam(fit_gam_vmh) # ~1min to predict
  saveRDS(pred_vmh, file = paste0("pred_gam_vmh_sample",sample,".rds") )
  
}

pred_vmh_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  pred_vmh  =readRDS( paste0("pred_gam_vmh_sample",sample,".rds") )
  pred_vmh_c[which_sample] = pred_vmh
}
set.seed(111)
domain_if_VMH = kmeans(pred_vmh_c, 2)$cluster
colData(spe_joint)$domain_if_VMH = domain_if_VMH
saveRDS(domain_if_VMH, file = paste0("domain_if_VMH.rds") )



domain_if_VMH_v2= array(dim=ncol(spe_joint))
for(sample in sample_names){
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  pred_vmh  =readRDS( paste0("pred_gam_vmh_sample",sample,".rds") )
  set.seed(111)
  domain_if_VMH_sample_tmp =  kmeans(pred_vmh, 2)$cluster
  domain_if_VMH_v2[which_sample] = domain_if_VMH_sample_tmp
}
# domain_if_VMH_v2 = kmeans(pred_vmh_c_v2, 2)$cluster
colData(spe_joint)$domain_if_VMH_v2 = domain_if_VMH_v2
saveRDS(domain_if_VMH_v2, file = paste0("domain_if_VMH_v2.rds") )
######################################################################### 
## do the plots
######################################################################### 
plt <- function(se_, cnames = "domain_if_VMH",
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){

  
  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  plot_bank1 <- plotColData(se_[,grp1],
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank2 <- plotColData(se_[,grp2],
                           x = "x_location", y = "y_location",
                           point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  
  plot_bank_byClu1 <- plotColData(se_[,grp1],
                                 x = "x_location", y = "y_location",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu2 <- plotColData(se_[,grp2],
                                 x = "x_location", y = "y_location",
                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
    scale_shape_manual(values=seq(0,15))
  
  
  # Visualize UMAPs of the non-spatial and BANKSY embedding:
  rdnames <- reducedDimNames(se_)
  
  umap_1 <- plotReducedDim(se_[,grp1],
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 2)+   
    scale_shape_manual(values=seq(0,15))
  
  umap_2<- plotReducedDim(se_[,grp2],
                          dimred = grep(paste0("UMAP.*lam", lambda, "$"), rdnames, value = TRUE),
                          colour_by = cnames, shape_by = "sample_id"
  ) + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda))+ 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~shape_by, ncol = 2)+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)

  print(plot_bank_byClu1)
  print(plot_bank_byClu2)
  
  print(umap_1)
  print(umap_2)
  
}

se = spe_joint
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])


pdf("multisample_banksy_2runs_13samples_domain_use_ARC_markerGenes.pdf", width = 60, height = 100/2)
plt(se[,], cnames = "domain_if_ARC")
dev.off()


pdf("multisample_banksy_2runs_13samples_domain_use_ARC_v2_markerGenes.pdf", width = 60, height = 100/2)
plt(se[,], cnames = "domain_if_ARC_v2")
dev.off()


pdf("multisample_banksy_2runs_13samples_domain_use_VMH_markerGenes.pdf", width = 60, height = 100/2)
plt(se[,], cnames = "domain_if_VMH")
dev.off()

pdf("multisample_banksy_2runs_13samples_domain_use_VMH_v2_markerGenes.pdf", width = 60, height = 100/2)
plt(se[,], cnames = "domain_if_VMH_v2")
dev.off()


