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

###########
# sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
###########
###########
sample_names_old =  c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")
sample_names_new =   c("br1225a","br1225b","br8741c",
                       "br8741d","br5459a","br5459b",
                       "br8667c")
sample_names= c(sample_names_old, sample_names_new)

########### load multi-sample Banksy result
spe_joint  = readRDS(paste0("spe_joint_banksy_lambda1_13samples_k_geom12.rds"))

# spe_joint_2runs_banksy_lambda0.rds
# se = spe_joint
# unique(colData(se)[,"sample_id"])
# colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
# unique(colData(se)[,"sample_id"] )


######################################################################### 
## find domains
######################################################################### 
# library(mgcv)
#### load results of smoothed domain after 2 rounds of knn smoothing
x1 = colData(spe_joint)$x_position
x2 = colData(spe_joint)$y_position

# knn_pred_ARC_c = array(dim=length(x1))
# for(sample in sample_names){
#   print(sample)
#   which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
#   knn_pred = readRDS(paste0("knn_predARC_",sample,"_k200_2ndRun.rds"))
#   # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,".rds") )
#   knn_pred_ARC_c[which_sample] = knn_pred
# }

# pred_ARC = (knn_pred_ARC_c>cutoff)

######################################################################### 
## find knn for each of within-domain cells
df_ = data.frame(x1=x1,x2=x2 )
library(FNN)

cutoff = 0.3

for(sample in sample_names){
  print(sample)
  
  # load results of smoothed domain after 2 rounds of knn smoothing
  knn_pred = readRDS(paste0("knn_predARC_",sample,"_k500_2ndRun.rds"))
  pred_ARC_which = which(knn_pred>cutoff)
  
  # find locations of cells within this sample and get knn
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  knn_ = get.knn(data=data.matrix(df_[which_sample, c("x1","x2")]), 
                 k=50)
  
  vec_cells_knn = unique(as.integer(knn_$nn.index[pred_ARC_which,]))
  
  vec_cells_knn_outsideDomain = vec_cells_knn[!vec_cells_knn  %in% pred_ARC_which]
  vec_cells_insideDomain = pred_ARC_which
  vec_cells_others = c(1:length(which_sample))[-c(vec_cells_knn_outsideDomain,pred_ARC_which)]
  saveRDS(vec_cells_knn_outsideDomain, file = paste0("vec_cells_knn_outsideARC_",sample,"_k500_2ndRun.rds") )
  saveRDS(vec_cells_insideDomain, file = paste0("vec_cells_insideARC_",sample,"_k500_2ndRun.rds") )
  saveRDS(vec_cells_others, file = paste0("vec_cells_others_ARC_",sample,"_k500_2ndRun.rds") )
  
  # knn_pred = unlist(lapply(1:length(which_sample), function(xx){
  # mean(knn_pred_ARC_c[which_sample[knn_$nn.index[xx,]]])
  # }))
  # saveRDS(knn_pred, file = paste0("knn_predARC_",sample,"_k500_2ndRun.rds") )
}



######################################################################### 
## do the plots
######################################################################### 
# get marker gene's expression
markerGene = "GLRA2"
vec_allSamples_exp_knn_outsideDomain = c()
vec_allSamples_exp_insideDomain  = c()
vec_allSamples_exp_others = c()

vec_allSamples_which_knn_outsideDomain = c()
vec_allSamples_which_insideDomain  = c()
vec_allSamples_which_others = c()

for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  
  vec_cells_knn_outsideDomain = readRDS(paste0("vec_cells_knn_outsideARC_",sample,"_k500_2ndRun.rds"))
  vec_cells_insideDomain = readRDS(paste0("vec_cells_insideARC_",sample,"_k500_2ndRun.rds"))
  vec_cells_others = readRDS(paste0("vec_cells_others_ARC_",sample,"_k500_2ndRun.rds"))
  
  exp_markergene_withinSample =  as.numeric((assays(spe_joint)$normcounts[markerGene,which_sample]))
  
  vec_allSamples_exp_knn_outsideDomain = c(vec_allSamples_exp_knn_outsideDomain, 
                                           exp_markergene_withinSample [vec_cells_knn_outsideDomain])
  
  vec_allSamples_exp_insideDomain = c(vec_allSamples_exp_insideDomain, exp_markergene_withinSample[vec_cells_insideDomain])
  vec_allSamples_exp_others = c(vec_allSamples_exp_others, exp_markergene_withinSample[vec_cells_others])
  
  
  vec_allSamples_which_knn_outsideDomain = c(vec_allSamples_which_knn_outsideDomain, which_sample[vec_cells_knn_outsideDomain] )
  vec_allSamples_which_insideDomain = c(vec_allSamples_which_insideDomain, which_sample[vec_cells_insideDomain] )
  vec_allSamples_which_others = c(vec_allSamples_which_others, which_sample[vec_cells_others] )
  
}


# do the plots - boxplt
df_ = data.frame(expr = c(vec_allSamples_exp_knn_outsideDomain,vec_allSamples_exp_others, vec_allSamples_exp_insideDomain ),
                 group = c(rep("near_domain", length(vec_allSamples_exp_knn_outsideDomain) ),
                           rep("far_from_domain", length(vec_allSamples_exp_others) ),
                           rep("within_domain", length( vec_allSamples_exp_insideDomain) )))
# "far_from_domain","within_domain"))
df_$group = factor(df_$group, levels = c("within_domain","near_domain","far_from_domain"))

summary(vec_allSamples_exp_insideDomain)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00   13.53   26.38   37.88  441.18 

pdf("validation_ARC_domain_detection.pdf", width = 15)
p = ggplot(df_, aes(x=group, y=expr, fill=group)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2,outlier.shape = NA) +
  # scale_fill_viridis(discrete = TRUE) +
  # theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("ARC") +
  xlab("")+ylim(c(0,40))
p
dev.off()


######################################################################### 
## do the plots - spatial plot
######################################################################### 
######################################################################### 
plt <- function(se_, cnames, 
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  
  colData(se_) <- cbind(colData(se_), spatialCoords(se_))
  
  plot_bank1 <- plotColData(se_[,grp1],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    # ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  plot_bank2 <- plotColData(se_[,grp2],
                            x = "x_location", y = "y_location",
                            point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  ) + coord_equal() + 
    # ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
    theme(plot.title = element_text(size=33)) + 
    facet_wrap(~ shape_by, ncol = 3)+   
    scale_shape_manual(values=seq(0,15))
  
  # 
  # plot_bank_byClu1 <- plotColData(se_[,grp1],
  #                                 x = "x_location", y = "y_location",
  #                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  # ) + coord_equal()+ 
  #   # ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
  #   theme(plot.title = element_text(size=33)) + 
  #   facet_wrap(~ colour_by + shape_by, ncol = length(sample_names_old))+   
  #   scale_shape_manual(values=seq(0,15))
  # 
  # plot_bank_byClu2 <- plotColData(se_[,grp2],
  #                                 x = "x_location", y = "y_location",
  #                                 point_size = 0.6, colour_by = cnames, shape_by = "sample_id"
  # ) + coord_equal()+ 
  #   # ggtitle(paste0("multi-sample Banksy, lambda = ", lambda,", cutoff = ",cutoff_)) + 
  #   theme(plot.title = element_text(size=33)) + 
  #   facet_wrap(~ colour_by + shape_by, ncol =  length(sample_names_new))+   
  #   scale_shape_manual(values=seq(0,15))
  
  print(plot_bank1)  
  print(plot_bank2)
  
  # print(plot_bank_byClu1)
  # print(plot_bank_byClu2)
  
}




se = spe_joint
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])


# vec_allSamples_which_knn_outsideDomain = c(vec_allSamples_which_knn_outsideDomain, which_sample[vec_cells_knn_outsideDomain] )
# vec_allSamples_which_insideDomain = c(vec_allSamples_which_insideDomain, which_sample[vec_cells_insideDomain] )
# vec_allSamples_which_others = c(vec_allSamples_which_others, which_sample[vec_cells_others] )
# df_$group = factor(df_$group, levels = c("within_domain","near_domain","far_from_domain"))

group_  = array(dim=ncol(se))
group_[vec_allSamples_which_knn_outsideDomain] = "near_domain"
group_[vec_allSamples_which_others] = "far_from_domain"
group_[vec_allSamples_which_insideDomain] = "within_domain"

colData(se)$group = group_

colData(se)$group= factor(colData(se)$group, levels = c("within_domain","near_domain","far_from_domain"))

pdf("validation_ARC_domain_detection_spatial.pdf", width = 60, height = 100/2)
plt(se[,], cnames = "group")
dev.off()



