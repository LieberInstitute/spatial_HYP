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



######################################################################### 
## plot 1st sample, for presentation -- cluster 33 and 25
######################################################################### 
coul <- brewer.pal(9, "Set1")
se = spe_joint

knn_pred_VMH_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predVMH_",sample,"_k200_2ndRun.rds"))
  knn_pred_VMH_c[which_sample] = knn_pred
}
knn_pred_VMH_c__ = knn_pred_VMH_c>0.2
names(knn_pred_VMH_c__)  =colData(spe_joint)$cell_id

# colData(se)$pred_VMH = knn_pred_VMH_c__
colData(se)$VMH = "others"
colData(se)$VMH[which(knn_pred_VMH_c__)] = "VMH"
colData(se)$VMH = factor(colData(se)$VMH, levels= c("others","VMH"))



######################################################################### 
## plot all samples, for presentation -- cluster 32 only, split by sex
######################################################################### 


# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
cellClu = as.character(colData(se)[,cnames])
cellClu[cellClu=="32"] = "cluster 32"
cellClu[cellClu!="cluster 32"] = "others"
cellClu[cellClu!="cluster 32" & knn_pred_VMH_c__] = "VMH"
cellClu = factor(cellClu, levels= c("others","VMH","cluster 32"))
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

pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu32_VMH_allSamples.pdf"), width = 20*7, height = 20*2)
plot_grid(plotlist = plot_bank, ncol = 7, byrow = TRUE)
dev.off()



######################################################################### 
## plot all samples, for presentation -- cluster 18 only, split by sex
######################################################################### 


# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)][1]
cellClu = as.character(colData(se)[,cnames])
cellClu[cellClu=="18"] = "cluster 18"
cellClu[cellClu!="cluster 18"] = "others"
cellClu[cellClu!="cluster 18" & knn_pred_VMH_c__] = "VMH"
cellClu = factor(cellClu, levels= c("others","VMH","cluster 18"))
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

pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu18_VMH_allSamples.pdf"), width = 20*7, height = 20*2)
plot_grid(plotlist = plot_bank, ncol = 7, byrow = TRUE)
dev.off()


######################################################################### 
## plot all samples, for presentation -- cluster 23 only, split by sex
######################################################################### 


# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)][1]
cellClu = as.character(colData(se)[,cnames])
cellClu[cellClu=="23"] = "cluster 23"
cellClu[cellClu!="cluster 23"] = "others"
cellClu[cellClu!="cluster 23" & knn_pred_VMH_c__] = "VMH"
cellClu = factor(cellClu, levels= c("others","VMH","cluster 23"))
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

pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu23_VMH_allSamples.pdf"), width = 20*7, height = 20*2)
plot_grid(plotlist = plot_bank, ncol = 7, byrow = TRUE)
dev.off()

######################################################################### 
## plot all samples, for presentation -- cluster 34 only, split by sex
######################################################################### 


# se_ = se[,colData(se)[,"sample_id"]=="br1735a"& colData(se)[,"clust_M0_lam0_k50_res2"] == "33"]
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)][1]
cellClu = as.character(colData(se)[,cnames])
cellClu[cellClu=="34"] = "cluster 34"
cellClu[cellClu!="cluster 34"] = "others"
cellClu[cellClu!="cluster 34" & knn_pred_VMH_c__] = "VMH"
cellClu = factor(cellClu, levels= c("others","VMH","cluster 34"))
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

pdf(paste0("multisample_2runs_banksy_lambda0_res", 0.6, "_for_presentation_clu34_VMH_allSamples.pdf"), width = 20*7, height = 20*2)
plot_grid(plotlist = plot_bank, ncol = 7, byrow = TRUE)
dev.off()
