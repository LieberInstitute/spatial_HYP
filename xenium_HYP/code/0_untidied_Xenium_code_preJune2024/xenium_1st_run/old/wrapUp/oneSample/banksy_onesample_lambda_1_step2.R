# ml conda_R/4.3
# R

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"

dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/wrapup/one_sample_banksy/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/wrapup/one_sample_banksy/"

dir.create("/dcs04/hansen/data/ywang/xenium/wrapup/")
dir.create("/dcs04/hansen/data/ywang/xenium/wrapup/one_sample_banksy")

setwd(dir_out)

###########
library(Seurat)
library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)

library(scater)
library(cowplot)
library(ggplot2)


######################################################################### 

xens2 = readRDS(paste0(dir_data_processed,
                       "xens2.rds")) # spe obj, not too much QC filteration

######################################################################### 
sampleID_unique = as.character(unique(xens2$sample_id))


ksample=1

# for(ksample in 1){
  sampleID_tmp = sampleID_unique[ksample]
  
  se = readRDS(paste0(dir_data_processed_banksy, 
                      "se_sample",sampleID_tmp,
                      "_preprosssed_lambda1_UMAP.rds"))
  
  
  
  ######################3
  for(res in c(0.5, 1, 1.2, 2, 4, 5)){
    
    print("res")
    print(res)
    se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = res)
    
    # Different clustering runs can be relabeled to minimise their differences with connectClusters:
    # se  <- Banksy::connectClusters(se)
    #> clust_M1_lam0.2_k50_res1.2 --> clust_M1_lam0_k50_res1.2
    #> 
    # Visualise the clustering output for non-spatial clustering (lambda=0) and BANKSY clustering (lambda=0.2).
    # cnames <- colnames(colData(se))
    # cnames <- cnames[grep("^clust", cnames)]
    colData(se) <- cbind(colData(se), spatialCoords(se))
    
    
    colData_ = colData(se)
    saveRDS(colData_, file = paste0(dir_data_processed_banksy, 
                                    "colData__sample",sampleID_tmp,"_preprosssed_res",
                                    res,"_lambda1.rds"))
    # saveRDS(se, file = paste0(dir_data_processed_banksy, 
                              # "se_sample",sampleID_tmp,"_preprosssed_res",
                              # res,"_lambda1.rds"))
    
  }
 
# }


# saveRDS(sampleID_unique, file=paste0(dir_data_processed_banksy, "sampleID_unique.rds"))

######################################################################### 
### plot the results
######################################################################### 
# plist = list()
# sampleID_unique  =readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
# 
# 
# for(ksample in 1){
#   for(res in c(1, 1.2, 2, 4, 5)){
#     
#     pdf(paste0("banksy_sample1_res",res,"_lambda1.pdf"), 
#         width = 10, height = 10)
#     
#     # saveRDS(se, file = paste0(dir_data_processed_banksy, 
#                               # "se_sample",sampleID_tmp,"_preprosssed_res",
#                               # res,".rds"))
#     # se = readRDS(paste0(dir_data_processed_banksy, 
#                         # "colData__sample",sampleID_tmp,"_preprosssed_res",
#                         # res,".rds"))
#     colData_ = readRDS(paste0(dir_data_processed_banksy, 
#                               "colData__sample",sampleID_tmp,"_preprosssed_res",
#                               res,"_lambda1.rds"))
#     colData(se) = colData_
#     
#     print(ksample)
#     sampleID_tmp = sampleID_unique[ksample]
#     se = readRDS(paste0(dir_data_processed_banksy, "se_sample",sampleID_tmp,"_preprosssed.rds"))
#     cnames <- colnames(colData(se))
#     cnames <- cnames[grep("^clust", cnames)]
#     colData(se) <- cbind(colData(se), spatialCoords(se))
#     
#     # plot_nsp <- plotColData(se,
#     #                         x = "x_location", y = "y_location",
#     #                         point_size = 0.6, colour_by = cnames[1]
#     # ) + coord_equal() + ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
#     
#     plot_bank <- plotColData(se,
#                              x = "x_location", y = "y_location",
#                              point_size = 0.6, colour_by = cnames[2]
#     ) + coord_equal()+ ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
#     
#     # plist[[ksample*2-1]] = plot_nsp
#     # plist[[ksample*2]] = plot_bank
#     
#     # Visualize UMAPs of the non-spatial and BANKSY embedding:
#     
#     rdnames <- reducedDimNames(se)
#     
#     # umap_nsp <- plotReducedDim(se,
#     # dimred = grep("UMAP.*lam$", rdnames, value = TRUE),
#     # colour_by = cnames[1]
#     # ) + ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
#     umap_bank <- plotReducedDim(se,
#                                 dimred = grep("UMAP.*lam1$", rdnames, value = TRUE),
#                                 colour_by = cnames[2]
#     ) + ggtitle(sampleID_tmp)+ theme(plot.title = element_text(size=33))
#     
#     print(plot_bank)
#     print(plot_bank + facet_wrap(~colour_by))
#     print(umap_bank)
#     
#     # print(plot_grid(plot_nsp , plot_bank , ncol = 2))
#     # print(plot_grid(
#     #   # plot_nsp + facet_wrap(~colour_by),
#     #   plot_bank + facet_wrap(~colour_by),
#     #   ncol = 2
#     # ))
#     # 
#     # print(plot_grid(
#     #   umap_nsp,
#     #   umap_bank,
#     #   ncol = 2
#     # ))
#   }
#   
#   dev.off()
#   
# }
# 
# 
# 





