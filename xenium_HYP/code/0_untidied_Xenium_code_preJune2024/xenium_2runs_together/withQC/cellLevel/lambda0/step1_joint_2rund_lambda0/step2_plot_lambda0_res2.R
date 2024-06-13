######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 0
res = 2
###########
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
# dir.create(dir_data_processed_banksy)

dir_data_processed = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
# dir_out = "/dcs04/hansen/data/ywang/xenium_newData/out/"
# setwd(dir_out)
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"
# dir.create(dir_data_processed_banksy)
# dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"

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
library(RColorBrewer)
###########
# sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
###########
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names= c(sample_names_old, sample_names_new)
sex_perSampleID = c("Male", "Female", "Female", "Male", "Male", "Female",
                    c("Male","Male","Female","Female","Male","Male","Female"))
names(sex_perSampleID) = sample_names
########### load multi-sample Banksy result
spe_joint = readRDS(paste0("spe_joint_2runs_banksy_lambda", lambda, "_res", res, ".rds"))
# spe_joint_2runs_banksy_lambda0.rds
se = spe_joint
unique(colData(se)[,"sample_id"])
colData(se)[,"sample_id"] = as.character(colData(se)[,"sample_id"])
unique(colData(se)[,"sample_id"] )

######################################################################### 
## do the plots
######################################################################### 
plt <- function(se_, 
                grp1 = colData(se)[,"sample_id"]%in%sample_names_old, 
                grp2 = colData(se)[,"sample_id"]%in%sample_names_new){
  
  cnames <- colnames(colData(se_))
  
  cnames <- cnames[grep("^clust", cnames)]
  print(cnames)
  
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
pdf(paste0("multisample_banksy_2runs_lambda_0_res",res,".pdf"), width = 60, height = 500/2)
plt(se[,])
dev.off()




######################################################################### 
## plot 1st sample, for presentation -- cluster 33 and 25
######################################################################### 
coul <- brewer.pal(9, "Set1")
se = spe_joint

knn_pred_arc_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predARC_",sample,"_k500_2ndRun.rds"))
  knn_pred_arc_c[which_sample] = knn_pred
}
knn_pred_arc_c__ = knn_pred_arc_c>0.3
names(knn_pred_arc_c__)  =colData(spe_joint)$cell_id

# colData(se)$pred_arc = knn_pred_arc_c__
colData(se)$ARC = "others"
colData(se)$ARC[which(knn_pred_arc_c__)] = "ARC"
colData(se)$ARC = factor(colData(se)$ARC, levels= c("others","ARC"))

###
se_ = se[,colData(se)[,"sample_id"]=="br1735a"]
knn_pred_arc_c__sampletmp = knn_pred_arc_c__[colData(se)[,"sample_id"]=="br1735a"]

xlim_ = range(se_$x_position)
ylim_ = range(se_$y_position)

# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se_))
cnames <- cnames[grep("^clust", cnames)]
cellClu = as.character(colData(se_)[,cnames])
cellClu[cellClu=="33"] = "cluster 33"
cellClu[cellClu!="cluster 33"] = "others"
cellClu[cellClu!="cluster 33" & knn_pred_arc_c__sampletmp] = "ARC"
cellClu = factor(cellClu, levels= c("others","ARC","cluster 33"))
# print(cnames)
colData(se_)$cluster = cellClu
# colData(se_) <- cbind(colData(se_), spatialCoords(se_))

# se = spe_joint

pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu33.pdf"), width = 20, height = 20)

# plot_bank1 <- plotColData(se_[,],
#                           x = "x_location", y = "y_location",
#                           point_size = 0.6, colour_by = "cluster", shape_by = "sample_id"
# ) + coord_equal() + xlim(xlim_)+ ylim(ylim_)+
#   ggtitle(paste0("multi-sample Banksy, lambda = ", lambda)) +
#   theme(plot.title = element_text(size=33),
#         axis.text = element_text(size=22),
#         axis.title=element_text(size=55)) +
#   # facet_wrap(~ shape_by, ncol = 3)+
#   scale_shape_manual(values=seq(0,15))
# plot_bank1
getPalette = colorRampPalette(coul)
palette = getPalette(nclu)[c(3,3,1)]

# nclu=3
# e41a1c
# pallet_ = c("grey","#377eb8","#4daf4a")

pallet_ = c("grey","#87b9e1","#e41a1c")
colData(se_)$y_position = - colData(se_)$y_position

plotSpots(se_[,], y_coord = "y_position", x_coord = "x_position",
          annotate = "cluster",

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


pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu25.pdf"), width = 20, height = 20)
# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se_))
cnames <- cnames[grep("^clust", cnames)][1]
cellClu = as.character(colData(se_)[,cnames])
cellClu[cellClu=="25"] = "cluster 25"
cellClu[cellClu!="cluster 25"] = "others"
cellClu[cellClu!="cluster 25" & knn_pred_arc_c__sampletmp] = "ARC"
cellClu = factor(cellClu, levels= c("others","ARC","cluster 25"))
# print(cnames)
colData(se_)$cluster = cellClu
getPalette = colorRampPalette(coul)
palette = getPalette(nclu)[c(3,3,1)]

# nclu=3
# e41a1c
# pallet_ = c("grey","#377eb8","#4daf4a")

pallet_ = c("grey","#87b9e1","#e41a1c")
colData(se_)$y_position = - colData(se_)$y_position

plotSpots(se_[,], y_coord = "y_position", x_coord = "x_position",
          annotate = "cluster",
          
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

# df <- cbind.data.frame(colData(se), spatialCoords(se))
# df <- cbind(colData(se), spatialCoords(se))

# 


######################################################################### 
## plot all samples, for presentation -- cluster 33 only, split by sex
######################################################################### 


# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
cellClu = as.character(colData(se)[,cnames])
cellClu[cellClu=="33"] = "cluster 33"
cellClu[cellClu!="cluster 33"] = "others"
cellClu[cellClu!="cluster 33" & knn_pred_arc_c_] = "ARC"
cellClu = factor(cellClu, levels= c("others","ARC","cluster 33"))
# print(cnames)
colData(se)$cluster = cellClu



# getPalette = colorRampPalette(coul)
# palette = getPalette(nclu)[c(3,3,1)]
pallet_ = c("grey","#87b9e1","#e41a1c")
colData(se_)$y_position = - colData(se_)$y_position
samples_kp = sample_names[sample_names!="br8667c"]
samples_kp = samples_kp[order(sex_perSampleID[samples_kp], decreasing = T)]

plot_bank <- lapply(samples_kp, function(x) {
  se_ = se[,colData(se)[,"sample_id"] == x]
  colData(se_)$y_position = - colData(se_)$y_position
  
  plotSpots(se_[,], y_coord = "y_position", x_coord = "x_position",
            annotate = "cluster",
            # annotate = cnames[1],
            size = 0.8,
            palette = pallet_,
            in_tissue = NULL) +
    theme(plot.title = element_text(size=33),legend.position = "none",
          axis.text = element_text(size=22),
          axis.title=element_text(size=55)) +   # labs(title = paste0(x, ",BANKSY clusters"))+
    theme(plot.title = element_text(size=33)) +
    coord_equal()

})

pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu33_allSamples.pdf"), width = 20*7, height = 20*2)
plot_grid(plotlist = plot_bank, ncol = 7, byrow = TRUE)
dev.off()



######################################################################### 
## plot all samples, for presentation -- cluster 33 only, split by sex
######################################################################### 


# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
cellClu = as.character(colData(se)[,cnames])
cellClu[cellClu=="25"] = "cluster 25"
cellClu[cellClu!="cluster 25"] = "others"
cellClu[cellClu!="cluster 25" & knn_pred_arc_c__] = "ARC"
cellClu = factor(cellClu, levels= c("others","ARC","cluster 25"))
# print(cnames)
colData(se)$cluster = cellClu


library(ggspavis)
# getPalette = colorRampPalette(coul)
# palette = getPalette(nclu)[c(3,3,1)]
pallet_ = c("grey","#87b9e1","#e41a1c")
colData(se_)$y_position = - colData(se_)$y_position
samples_kp = sample_names[sample_names!="br8667c"]
samples_kp = samples_kp[order(sex_perSampleID[samples_kp], decreasing = T)]

plot_bank <- lapply(samples_kp, function(x) {
  se_ = se[,colData(se)[,"sample_id"] == x]
  colData(se_)$y_position = - colData(se_)$y_position
  
  plotSpots(se_[,], y_coord = "y_position", x_coord = "x_position",
            annotate = "cluster",
            # annotate = cnames[1],
            size = 0.8,
            palette = pallet_,
            in_tissue = NULL) +
    theme(plot.title = element_text(size=33),legend.position = "none",
          axis.text = element_text(size=22),
          axis.title=element_text(size=55)) +   # labs(title = paste0(x, ",BANKSY clusters"))+
    theme(plot.title = element_text(size=33)) +
    coord_equal()
  
})

pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu25_allSamples.pdf"), width = 20*7, height = 20*2)
plot_grid(plotlist = plot_bank, ncol = 7, byrow = TRUE)
dev.off()

