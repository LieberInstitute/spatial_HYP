######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
# dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
# setwd(dir_out)
###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
# dir.create(dir_data_processed_banksy)
dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/")
# dir.create("/dcs04/hansen/data/ywang/xenium_joint_runs/out/")

setwd(dir_out)

###########
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Banksy)
aname <- "normcounts"
library(scater)

###########
sample_names_old = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))
sample_names_new = readRDS(paste0(dir_data_processed_banksy_new, "sampleID_unique.rds"))

########### combine 2 runs, using ARC cells only
list_subdomains_arc_ARConly_combine2Runs = list()

list_clusID_subdomains_arc_ARConly_combine2Runs = list()
list_clusID_subdomains_arc_ARConly_combine2Runs[["subdomain_1"]] = list(res1 = c(9,11), res2= c(16, 18))
list_clusID_subdomains_arc_ARConly_combine2Runs[["subdomain_2"]] = list(res0.6 = 5 , res1 = 10, res2= c( 7))
# list_clusID_subdomains_arc_ARConly_combine2Runs[["subdomain_2"]] = list(res0.6 = 5 , res1 = 10, res2= c(20, 11, 7))
list_clusID_subdomains_arc_ARConly_combine2Runs[["subdomain_3"]] =  list(res0.6 = 2 , res1 = 2, res2= c(17, 12))
# list_clusID_subdomains_arc_ARConly_combine2Runs[["subdomain_3"]] =  list(res0.6 = 2 , res1 = 2, res2= c(17, 12, 13))

# list_clusID_domain2_arc_ARConly_combine2Runs = list(res0.6 = 5 , res1 = 10, res2= c(20, 11, 7))
# list_clusID_domain3_arc_ARConly_combine2Runs = list(res0.6 = 2 , res1 = 2, res2= c(17, 12, 13))
for(subdomain in 1:3){
  list_tmp = list()
  
  for(res in c( 0.6, 1, 2)){
    print(res)
    colData_ = readRDS( paste0("colData_spe_joint_2runs_banksy_lambda0_ARConly_res",res,".rds"))
    # saveRDS(colData_, file = paste0("colData_spe_joint_2runs_banksy_lambda0_ARConly_res",res,".rds"))
    clus_tmp = colData_[,paste0("clust_M0_lam0_k50_res", res)]
    clus_withinSubdomain = list_clusID_subdomains_arc_ARConly_combine2Runs[[paste0("subdomain_", subdomain)]][[paste0("res",res)]]
    list_tmp[[paste0("res",res)]] = which(clus_tmp %in% clus_withinSubdomain)
    
  }
  
  list_subdomains_arc_ARConly_combine2Runs[[paste0("subdomain_", subdomain)]] = list_tmp
}
saveRDS(list_subdomains_arc_ARConly_combine2Runs, file = paste0("list_subdomains_arc_ARConly_combine2Runs.rds"))

# summarize table for this setting, comparing across resolutions
list_table_ = list()
for(subdomain in 1:3){
  print(subdomain)
  table_tmp = array(dim=c(3,3))
  colnames(table_tmp) = rownames(table_tmp) = paste0("res", c( 0.6, 1, 2))
  
  
  for(res1 in c( 0.6, 1, 2)){
    print(res1)
    
    which_cells_res1 =  list_subdomains_arc_ARConly_combine2Runs[[paste0("subdomain_", subdomain)]] [[paste0("res",res1)]]
    ncell1 = length(which_cells_res1)
    print(ncell1)
    
    for(res2 in c( 0.6, 1, 2)){
      which_cells_res2 =  list_subdomains_arc_ARConly_combine2Runs[[paste0("subdomain_", subdomain)]] [[paste0("res",res2)]]
      ncell2 = length(which_cells_res2)
      # print(ncell2)
      
      table_tmp[paste0("res", res1), paste0("res", res2)] = length(intersect(which_cells_res1,which_cells_res2))/mean(c(ncell1, ncell2))
    }
  }
  
  list_table_[[paste0("subdomain_", subdomain)]] = table_tmp
  
}

list_table_

subdomain=1
res1=1
res2=0.6

which_cells_res1 =  list_subdomains_arc_ARConly_combine2Runs[[paste0("subdomain_", subdomain)]] [[paste0("res",res1)]]
which_cells_res2 =  list_subdomains_arc_ARConly_combine2Runs[[paste0("subdomain_", subdomain)]] [[paste0("res",res2)]]


