######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R

###########
lambda <- 1
lambdaName="1"
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_data_processed_banksy = "/dcs04/hansen/data/ywang/xenium/data_processed_banksy/"
dir_out = "/dcs04/hansen/data/ywang/xenium/banksy/"
setwd(dir_out)

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
sample_names = readRDS(paste0(dir_data_processed_banksy, "sampleID_unique.rds"))

########### 
sampleID_unique <- c("br6197","br5993a","br5993b",
                     "br1735a","br1735b","br6588")
sex_perSampleID = c("M", "F", "F", "M", "M", "F")
names(sex_perSampleID) = sampleID_unique
# table(sex)
########### load multi-sample Banksy result
# spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,".rds"))

spe_joint  = readRDS(paste0("spe_joint_banksy_lambda",lambdaName,"_arcOnly.rds"))
# spe_joint_banksy_lambda1_res2.rds
# se = spe_joint
# Idents(se_seurat) = colData(se)[,paste0("clust_M0_lam", lambda, "_k50_res0.6")]

# head(colData(se))

colData(spe_joint)$clust_M0_lam1_k50_res0.6 = NULL
# colnames(colData(se))

# colData(se)[,which(colnames(colData(se))=="clust_M0_lam1_k50_res0.6")[1]] = NULL # rm the layer annotations using full data

###########
clu = colData(spe_joint)$clust_M0_lam1_k50_res0.6
#  use clu 1, clu 2+3, clu7, and else
clu_trans = as.integer(clu)
clu_trans[clu_trans==3]=2
clu_trans[!clu_trans%in%c(1,2,3,7)]=0

table(clu_trans)
# clu_trans
# 0     1     2     7 
# 26637 10341 15932  4799 

################### smooth within each sample and plot fitted clu1
# library(fields)
library(mgcv)
library(ggplot2)

get_smoothed_sp <- function(clu_bin, pos){
  
  df_=data.frame(x1=pos$x1,x2=pos$x2,clu_bin=clu_bin)
  fit_ <- gam(clu_bin~s(x1,x2),data=df_)
  
  pred=predict(fit_,
               df_)
  list(fit_=fit_,
       pred=pred )
}

ksample = 1
clu_tmp = 0

clus_unique = c(0, 1, 2, 7)
# df_smoothed_clu = array(dim=c(6,4))

for(ksample in 1:6){
  sample_tmp = sample_names[ksample]
  which_sample_tmp = which(colData(spe_joint)$sample_id==sample_tmp)
              
  pos = (colData(spe_joint)[which_sample_tmp,c("x_position","y_position")])
    colnames(pos)=c("x1","x2")
  for(clu_tmp in clus_unique){
      clu_tmp_bin=clu_trans
      clu_tmp_bin[clu_trans==clu_tmp]=1
      clu_tmp_bin[clu_trans!=clu_tmp]=0
      
      clu_tmp_smoothed_out = get_smoothed_sp(clu_tmp_bin[which_sample_tmp],
                                                  pos )
      df_tmp=data.frame(x1=pos$x1,
                   x2=pos$x2,
                   clu_tmp_smoothed=clu_tmp_smoothed_out$pred,
                   sample = sample_tmp)
                   # clu_tmp = clu_tmp)
      names(df_tmp)[names(df_tmp)=="clu_tmp_smoothed"] = paste0("cluster ",clu_tmp)
      
      if(clu_tmp == clus_unique[1]){
        df_sampleTmp = df_tmp
      }else{
        df_sampleTmp = cbind(df_sampleTmp, df_tmp[,3])
        colnames(df_sampleTmp)[(colnames(df_sampleTmp)=="df_tmp[, 3]")] = paste0("cluster ",clu_tmp)
      }
  

  }
  if(ksample == 1 ){
      df_ = df_sampleTmp
  }else{
      df_ = rbind(df_, df_sampleTmp) 
  }   
  # ggplot(df_tmp, aes(x1,x2,color=clu_tmp_smoothed)) +
    # geom_point()+
    # facet_wrap(~sample, ncol =3)
}
saveRDS(df_, file="df_smootedClu_arcOnly_.rds")

# saveRDS(df_, file="df_smootedClu_arcOnly.rds")


head(df_)

### reassign domain
set.seed(111)
domain_smoothed = kmeans(df_[,c("cluster 0", "cluster 1","cluster 2","cluster 7")],
                         4)$cluster

df_$domain_smoothed = as.character(domain_smoothed)
ggplot(df_, aes(x1,x2,color=domain_smoothed)) +
  geom_point()+
  facet_wrap(~sample, ncol =3)

# sample_tmp = sample_names[ksample]

for(kclu in 1:4){
  clu_tmp = clus_unique[kclu]
  df_clu_tmp = df_[,]
  df_clu_tmp$cluster = clu_tmp
  df_clu_tmp$fitted_prop = df_[,paste0("cluster ", clu_tmp)]
    
  if(kclu ==1){
    df_v2 = df_clu_tmp[, c("fitted_prop", "sample", "domain_smoothed","x1","x2","cluster")]
  }else{
    df_v2 = rbind(df_v2, df_clu_tmp[, c("fitted_prop", "sample", "domain_smoothed","x1","x2","cluster")])
  }
}
saveRDS(df_v2, file="df_v2_smootedClu_arcOnly_.rds")

df_v2_domainSmoothedOnly = df_v2[df_v2$cluster == clus_unique[1],]
df_v2_domainSmoothedOnly$cluster=NULL

saveRDS(df_v2_domainSmoothedOnly, file="df_v2_domainSmoothedOnly_arcOnly_.rds")

##########
library(ggplot2)=
## plot the new-defined sub-domains
df_v2_domainSmoothedOnly  = readRDS("df_v2_domainSmoothedOnly_arcOnly_.rds")
df_v2_domainSmoothedOnly$sex = sex_perSampleID[df_v2_domainSmoothedOnly$sample]
df_v2_domainSmoothedOnly$sex[df_v2_domainSmoothedOnly$sex=="M"] = "Male"
df_v2_domainSmoothedOnly$sex[df_v2_domainSmoothedOnly$sex=="F"] = "Female"
saveRDS(df_v2_domainSmoothedOnly, file="df_v2_domainSmoothedOnly_arcOnly_.rds")


png("smootedClu_arcOnly.png", height = 5, width = 10, res=400, units = "in")
ggplot(df_v2_domainSmoothedOnly, aes(x1,x2,color=domain_smoothed)) +
  geom_point(size=0.002)+
  facet_wrap(~sex+sample, ncol =3)
dev.off()

# scale_color_manual()

# plot(pos[,1],pos[,2],col=clu1_smoothed_out$pred)
# b <- getViz(clu1_smoothed_out$fit_)
# plot(sm(b, 1)) + l_fitRaster() + l_fitContour() + l_points()

