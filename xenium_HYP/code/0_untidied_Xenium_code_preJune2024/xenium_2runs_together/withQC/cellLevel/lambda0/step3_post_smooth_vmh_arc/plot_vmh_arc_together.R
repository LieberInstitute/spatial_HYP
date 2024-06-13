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
aname <- "normcounts"
library(ggspavis)

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
knn_pred_vmh_c_ = readRDS("knn_VMH_twoRounds.rds")
knn_pred_arc_c_ = readRDS("knn_ARC_twoRounds.rds")

######################################################################### 
## do the plots
######################################################################### 

######################################################################### 
## plot 1st sample, for presentation
######################################################################### 
pdf(paste0("knn_vmh_arc_br1735a.pdf"), width = 20, height = 20)
knn_pred_vmh_c__sampletmp = knn_pred_vmh_c_[colData(se)[,"sample_id"]=="br1735a"]
knn_pred_arc_c__sampletmp = knn_pred_arc_c_[colData(se)[,"sample_id"]=="br1735a"]
se_ = se[,colData(se)[,"sample_id"]=="br1735a"]
cnames = "domain"
colData(se_)$domain = "others"
colData(se_)$domain[knn_pred_vmh_c__sampletmp=="VMH"] = "VMH"
colData(se_)$domain[knn_pred_arc_c__sampletmp=="ARC"] = "ARC"
pallet_ = c("grey","#5DD959","#8C63CF")
names(pallet_) = c("others", "VMH", "ARC")
colData(se_)$y_position = - colData(se_)$y_position
colData(se_)$domain = factor(colData(se_)$domain, levels = c("ARC", "VMH", "others"))
plotSpots(se_[,], y_coord = "y_position", x_coord = "x_position",
          annotate = cnames,
          size = 0.8,
          palette = pallet_,
          in_tissue = NULL) +
  theme(plot.title = element_text(size=33),
        axis.text = element_text(size=22),
        axis.title=element_text(size=55)) +   # labs(title = paste0(x, ",BANKSY clusters"))+
  theme(plot.title = element_text(size=33)) +
  coord_equal()
dev.off()


pdf(paste0("knn_vmh_arc_br6197.pdf"), width = 20, height = 20)
knn_pred_vmh_c__sampletmp = knn_pred_vmh_c_[colData(se)[,"sample_id"]=="br6197"]
knn_pred_arc_c__sampletmp = knn_pred_arc_c_[colData(se)[,"sample_id"]=="br6197"]
se_ = se[,colData(se)[,"sample_id"]=="br6197"]
cnames = "domain"
colData(se_)$domain = "others"
colData(se_)$domain[knn_pred_vmh_c__sampletmp=="VMH"] = "VMH"
colData(se_)$domain[knn_pred_arc_c__sampletmp=="ARC"] = "ARC"
pallet_ = c("grey","#5DD959","#8C63CF")
names(pallet_) = c("others", "VMH", "ARC")
colData(se_)$y_position = - colData(se_)$y_position
colData(se_)$domain = factor(colData(se_)$domain, levels = c("ARC", "VMH", "others"))
plotSpots(se_[,], y_coord = "y_position", x_coord = "x_position",
          annotate = cnames,
          size = 0.8,
          palette = pallet_,
          in_tissue = NULL) +
  theme(plot.title = element_text(size=33),
        axis.text = element_text(size=22),
        axis.title=element_text(size=55)) +   # labs(title = paste0(x, ",BANKSY clusters"))+
  theme(plot.title = element_text(size=33)) +
  coord_equal()
dev.off()