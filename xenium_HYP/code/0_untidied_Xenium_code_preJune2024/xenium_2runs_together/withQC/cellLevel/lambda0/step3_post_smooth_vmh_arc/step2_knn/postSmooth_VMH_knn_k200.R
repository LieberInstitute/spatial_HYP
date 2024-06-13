######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 0

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
library(ggspavis)

aname <- "normcounts"
###########
# sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
###########
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_2runs_banksy_lambda0.rds"))

# spe_joint_2runs_banksy_lambda0.rds
se = spe_joint
unique(colData(se)[,"sample_id"])
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )


######################################################################### 
## find domains
######################################################################### 
# library(mgcv)
x1 = colData(spe_joint)$x_position
x2 = colData(spe_joint)$y_position

###########  smooth VMH and ARC
library(FNN)
# library("fastknn")
# library("caTools")

# ARC: 
# load qst run of knn results: k=50, cuoff = 0.1 
knn_pred_arc_c = array(dim=length(x1))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predARC_",sample,"_k200.rds"))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
  knn_pred_arc_c[which_sample] = knn_pred
}
knn_pred_arc = knn_pred_arc_c > 0.1

# run 2nd run of knn
df_ = data.frame(x1=x1,x2=x2, if_ARC = knn_pred_arc)

for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)

  knn_ = get.knn(data=data.matrix(df_[which_sample, c("x1","x2")]), 
                 k=50)
  
  knn_pred = unlist(lapply(1:length(which_sample), function(xx){
    mean(knn_pred_arc[which_sample[knn_$nn.index[xx,]]])
  }))
  
  saveRDS(knn_pred, file = paste0("knn_predARC_",sample,"_k200_2ndRun.rds") )
}


knn_pred_arc_c = array(dim=length(knn_pred_arc))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predARC_",sample,"_k200_2ndRun.rds"))
  knn_pred_arc_c[which_sample] = knn_pred
}



############################################################
#################### VMH
############################################################
# load result of 1st run knn on vmh: k=10, cutoff = 0.2
knn_pred_vmh_c = array(dim=length(x1))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predVMH_",sample,".rds"))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
  knn_pred_vmh_c[which_sample] = knn_pred
}
knn_pred_vmh = knn_pred_vmh_c > 0.2

# run the 2nd run of knn with k=10
df_ = data.frame(x1=x1,x2=x2, if_ARC = knn_pred_vmh_c)
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  knn_ = get.knn(data=data.matrix(df_[which_sample, c("x1","x2")]), 
                 k=200)
  
  knn_pred = unlist(lapply(1:length(which_sample), function(xx){
    mean(knn_pred_vmh[which_sample[knn_$nn.index[xx,]]])
  }))
  
  saveRDS(knn_pred, file = paste0("knn_predVMH_",sample,"_k200_2ndRun.rds") )
  
}


knn_pred_vmh_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predVMH_",sample,"_k200_2ndRun.rds"))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
  knn_pred_vmh_c[which_sample] = knn_pred
}

knn_pred_vmh_c_ = knn_pred_vmh_c>0.2
names(knn_pred_vmh_c_)  =colData(spe_joint)$cell_id
knn_pred_vmh_c_ = as.character(knn_pred_vmh_c_)
knn_pred_vmh_c_[knn_pred_vmh_c_=="FALSE"] = "not VMH"
knn_pred_vmh_c_[knn_pred_vmh_c_=="TRUE"] = "VMH"

# saveRDS(knn_pred_vmh_c_, file = paste0("knn_VMH_twoRounds.rds") )

######################################################################### 
## do the plots
######################################################################### 
######################################################################### 
plt <- function(se_, cnames, cutoff_,
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  
  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  plot_bank1 <- plotColData(se_[,grp1],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank2 <- plotColData(se_[,grp2],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  
  plot_bank_byClu1 <- plotColData(se_[,grp1],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank_byClu2 <- plotColData(se_[,grp2],
                                  x = "x_location", y = "y_location",
                                  point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal()+ 
    ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
    scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)
  
  print(plot_bank_byClu1)
  print(plot_bank_byClu2)
  
}




se = spe_joint
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])


pdf("multisample_banksy_2runs_13samples_domain_knn_pred_vmh_c_k200_2ndRun.pdf", width = 60, height = 100/2)
for(cutoff in c(0.05,0.1,0.2,0.3,0.4,0.5)){
  print(cutoff)
  colData(se)$pred_vmh = as.character(knn_pred_vmh_c>cutoff)
  
  plt(se[,], cnames = "pred_vmh",cutoff_=cutoff)
}
dev.off()

pdf("multisample_banksy_2runs_13samples_domain_knn_pred_vmh_c_k200_2ndRun_cutoff06_09.pdf", width = 60, height = 100/2)
for(cutoff in c(0.6, 0.7, 0.8,0.9)){
  print(cutoff)
  colData(se)$pred_vmh = as.character(knn_pred_vmh_c>cutoff)
  
  plt(se[,], cnames = "pred_vmh",cutoff_=cutoff)
}
dev.off()






######################################################################### 
## plot 1st sample, for presentation
######################################################################### 

pdf(paste0("knn_pred_vmh_for_presentation.pdf"), width = 20, height = 20)
cutoff = 0.2
# plot_grid(plotlist = plot_bank, ncol = 3, byrow = TRUE)
# plt()
print(cutoff)
# colData(se)$pred_vmh = (knn_pred_vmh_c>cutoff)
knn_pred_vmh_c__sampletmp = (knn_pred_vmh_c>cutoff)[colData(se)[,"sample_id"]=="br1735a"]

# plt(se[,], cnames = "pred_vmh",cutoff_=cutoff)
se_ = se[,colData(se)[,"sample_id"]=="br1735a"]
# cnames <- colnames(colData(se_))
# cnames <- cnames[grep("^clust", cnames)]
# print(cnames)
cnames = "VMH"
# colData(se_) <- cbind(colData(se_), spatialCoords(se_))
colData(se_)$VMH = "others"
colData(se_)$VMH[knn_pred_vmh_c__sampletmp] = "VMH"

# plotColData(se_[,],
#             x = "x_location", y = "y_location",
#             point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
# ) + coord_equal() + 
#   ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff)) + 
#   theme(plot.title = element_text(size=33),
#         axis.text = element_text(size=22),
#         axis.title=element_text(size=55)) + 
#   # facet_wrap(~ shape_by, ncol = 3)+   
#   scale_shape_manual(values=seq(0,15))
pallet_ = c("grey","#87b9e1")
colData(se_)$y_position = - colData(se_)$y_position

plotSpots(se_[,], y_coord = "y_position", x_coord = "x_position",
          annotate = cnames,
          
          # annotate = cnames[1],
          size = 0.8,
          # palette=coul,
          palette = pallet_,
          in_tissue = NULL) +
  theme(plot.title = element_text(size=33),
        axis.text = element_text(size=22),
        axis.title=element_text(size=55)) +   # labs(title = paste0(x, ",BANKSY clusters"))+
  theme(plot.title = element_text(size=33)) +
  coord_equal()
dev.off()
