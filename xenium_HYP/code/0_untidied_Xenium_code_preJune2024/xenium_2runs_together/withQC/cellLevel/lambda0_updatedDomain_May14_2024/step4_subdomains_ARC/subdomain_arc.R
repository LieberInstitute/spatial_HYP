######################################################################### 
# Running BANKSY
######################################################################### 
# ml conda_R/4.3
# R
###########
k=500
cutoff = 0.3
###########
lambdaName=0
lambda=0
res=2
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed_banksy =  "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream/"
###########
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_downstream_newRun/"

###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)

########### 
dir_data_processed_banksy_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed_banksy/"
dir_data_processed_new = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
sample_names_old  <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")

sample_names_new = c("br1225a","br1225b","br8741c",
                     "br8741d","br5459a","br5459b",
                     "br8667c")
sample_names = c(sample_names_old, sample_names_new)
########### load spe_joint

spe_joint = readRDS(paste0("spe_joint_2runs_banksy_lambda0_res", res, ".rds"))
########### load ARC domain identification result

knn_pred_arc_c = array(dim=ncol(spe_joint))
for(sample in sample_names){
  print(sample)
  which_sample = which(colData(spe_joint)[,"sample_id"] == sample)
  knn_pred = readRDS(paste0("knn_predARC_",sample,"_k500_2ndRun.rds"))
  knn_pred_arc_c[which_sample] = knn_pred
}

pred_arc = (knn_pred_arc_c > cutoff)

########## only keep the cells within ARC
spe_joint_arc= spe_joint[,pred_arc]

########### load cell type clustering result of ARC

cnames <- colnames(colData(spe_joint_arc))
cnames <- cnames[grep("^clust", cnames)]
print(cnames)

########## specify cell clusters with within-arc spatial variance
clu_sele = c("12","17","33","25")
# clu_sele = as.integer(clu_sele)
########## 
k_subdomain = 3

###########
clu = colData(spe_joint_arc)[,cnames]
clu = as.character(clu)
#  use clu 1, clu 2+3, clu7, and else
# clu = as.integer(clu)
# clu[clu==3]=2
clu[!clu%in%clu_sele]="0"
# clu = clu
# table(clu)
# clu
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

# ksample = 1
# clu_tmp = 0

# clu_sele = c(0, 1, 2, 7)
# df_smoothed_clu = array(dim=c(6,4))

for(ksample in 1:length(sample_names)){
  print("ksample")
  print(ksample)
  
  sample_tmp = sample_names[ksample]
  which_sample_tmp = which(colData(spe_joint_arc)$sample_id==sample_tmp)
  
  pos = (colData(spe_joint_arc)[which_sample_tmp,c("x_position","y_position")])
  colnames(pos)=c("x1","x2")
  for(clu_tmp in clu_sele){
    print("clu_tmp")
    print(clu_tmp)
    
    clu_tmp_bin=clu
    clu_tmp_bin[clu==clu_tmp]=1
    clu_tmp_bin[clu!=clu_tmp]=0
    clu_tmp_bin = as.integer(clu_tmp_bin)
    clu_tmp_smoothed_out = get_smoothed_sp(clu_tmp_bin[which_sample_tmp],
                                           pos )
    df_tmp=data.frame(x1=pos$x1,
                      x2=pos$x2,
                      clu_tmp_smoothed=clu_tmp_smoothed_out$pred,
                      sample = sample_tmp)
    # clu_tmp = clu_tmp)
    names(df_tmp)[names(df_tmp)=="clu_tmp_smoothed"] = paste0("cluster ",clu_tmp)
    
    if(clu_tmp == clu_sele[1]){
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
head(df_)

### reassign domain
set.seed(111)
domain_smoothed = kmeans(df_[,paste0("cluster ", clu_sele)],
                         3)$cluster

df_$domain_smoothed = as.character(domain_smoothed)
# ggplot(df_, aes(x1,x2,color=domain_smoothed)) +
#   geom_point()+
#   facet_wrap(~sample, ncol =3)

# sample_tmp = sample_names[ksample]

for(kclu in 1:length(clu_sele)){
  clu_tmp = clu_sele[kclu]
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


df_v2_domainSmoothedOnly = df_v2[df_v2$cluster == clu_sele[1],]
df_v2_domainSmoothedOnly$cluster=NULL

saveRDS(df_v2_domainSmoothedOnly, file="df_v2_domainSmoothedOnly_arcOnly_.rds")

##########
library(ggplot2)
## plot the new-defined sub-domains
df_v2_domainSmoothedOnly  = readRDS("df_v2_domainSmoothedOnly_arcOnly_.rds")
# df_v2_domainSmoothedOnly$sex = sex_perSampleID[df_v2_domainSmoothedOnly$sample]
# df_v2_domainSmoothedOnly$sex[df_v2_domainSmoothedOnly$sex=="M"] = "Male"
# df_v2_domainSmoothedOnly$sex[df_v2_domainSmoothedOnly$sex=="F"] = "Female"
# saveRDS(df_v2_domainSmoothedOnly, file="df_v2_domainSmoothedOnly_arcOnly_.rds")


png("smootedClu_arcOnly_k3.png", height = 15, width = 10, res=400, units = "in")
ggplot(df_v2_domainSmoothedOnly, aes(x1,x2,color=domain_smoothed)) +
  geom_point(size=0.002)+
  # facet_wrap(~sex+sample, ncol =3)
facet_wrap(~sample, ncol =3)
dev.off()




