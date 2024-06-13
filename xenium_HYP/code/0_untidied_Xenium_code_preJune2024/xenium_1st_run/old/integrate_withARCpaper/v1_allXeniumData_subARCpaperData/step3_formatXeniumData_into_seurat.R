# dir_Data_arcPaper = "/dcs04/hansen/data/ywang/xenium/data_paperARC/"
# 
dir_out= "//dcs04/hansen/data/ywang/xenium//integrate_paperARCscRNA/"
setwd(dir_out)
# 
# ########## load seurat object for scRNA data in ARC paper 
# obj_arcPaper = readRDS(paste0(dir_Data_arcPaper, "obj_arcPaper_sub.rds"))

########## load xenium data
lambdaName = 0
spe_joint  = readRDS(paste0("/dcs04/hansen/data/ywang/xenium/banksy/spe_joint_banksy_lambda",lambdaName,"_nucleus.rds"))
se = spe_joint

########## format xenium data
lambda=0
se <- scater::computeLibraryFactors(se)
se <- scater::logNormCounts(se)

se_seurat  <- as.Seurat(se) 
Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]
# se_seurat = FindAllMarkers(se_seurat)
se_seurat = NormalizeData(se_seurat)
se_seurat = FindVariableFeatures(se_seurat)
saveRDS(se_seurat, file="obj_xenium_nuclei.rds")