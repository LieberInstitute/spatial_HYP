# ml conda_R/4.3.x
# R
####################################
####################################
####################################

dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20240223__174743__022324_KMon120823/"
# dir_data = "/dcs04/hansen/data/ywang/xenium/raw_Data/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC_newRun/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out_QC_newRun/"
# dir.create(dir_data_processed)
# dir.create(dir_out)
# dir.create(paste0(dir_data_processed,"/cellLevel"))
# dir.create(paste0(dir_data_processed,"/nucleusLevel"))
# 
####################################
####################################
####################################
library(data.table)
library(Biostrings)
library(ggplot2)
library(gridExtra)
#library(SpatialExperiment)
library(spatialLIBD)
require(colorout)
library(default)
library(MoleculeExperiment)
library(magrittr)
library(scater)
library(scran)
library(BiocParallel)
library(job)
library(SpatialExperiment)
##################################################################################################
##################################################################################################
##################################################################################################
ColorOut()
options("styler.addins_style_transformer" = "biocthis::bioc_style()")

## change write.table to output a TSV without quotes or rownames by default
default(write.table) <- list(row.names=F,col.names=T,sep='\t',quote=F)

theme_set(theme_bw()+theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.y = element_text(size =16), plot.title = element_text(size = 20,hjust=0.5), strip.text = element_text(size=18), legend.text = element_text(size=10), legend.title = element_text(size=11,hjust=0.5)))



##################################################################################################
######## for exploration of the overall data qualtiy
##################################################################################################
# library(fread)
samples_filderName = list.files(paste0(dir_data))
for(sample_tmp in samples_filderName[]){
  print(sample_tmp)
  ## read transcripts.csv 
  # gunzip -k transcripts.csv.gz
  df_transcripts <- fread(paste0(dir_data,sample_tmp,
                                 # "/output-XETG00089__00011586__Region_3__20231019__161214/",
                                 "/transcripts.csv.gz"), header = T)
  ###### filter the transcripts
  # df_transcripts_QC20 = df_transcripts[df_transcripts$qv >= 20,]
  print(nrow(df_transcripts))
  print(sum(df_transcripts$qv >= 20))
  print(sum(df_transcripts$qv >= 20)/nrow(df_transcripts))
  
  ###### create a csv file with only filtered transcripts
  # dir.create(paste0(dir_data_processed,"/cellLevel/",sample_tmp,"_cellLevel_qv20"))
  
  # print("writing csv file")
  
  # fwrite(df_transcripts_QC20, file= paste0(dir_data_processed,"/cellLevel/",sample_tmp,"_cellLevel_qv20/", "transcripts.csv"))
  
  # ###### keep with-nuclei data
  # dir.create(paste0(dir_data_processed,"/nucleusLevel/",sample_tmp,"_nucleusLevel_qv20"))
  # 
  # df_transcripts_QC20_nuclear = df_transcripts_QC20[df_transcripts_QC20$overlaps_nucleus ==1  ,]
  # fwrite(df_transcripts_QC20_nuclear, file= paste0(dir_data_processed,"/nucleusLevel/",sample_tmp,"_nucleusLevel_qv20/", "transcripts.csv"))
  
}



# [1] "output-XETG00089__0011199__11199X_1225A_HYP__20240223__174831"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 52864949
# [1] 43720008
# [1] 0.8270132
# [1] "output-XETG00089__0011199__11199X_1225B_HYP__20240223__174831"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 36673547
# [1] 30204005
# [1] 0.823591
# [1] "output-XETG00089__0011199__11199X_8741C_HYP__20240223__174831"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 29207146
# [1] 23616715
# [1] 0.8085937
# [1] "output-XETG00089__0011199__11199X_8741D_HYP__20240223__174831"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 26150164
# [1] 21384922
# [1] 0.8177739
# [1] "output-XETG00089__0011236__11236X_5459A_HYP__20240223__174831"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 36342245
# [1] 31053908
# [1] 0.8544851
# [1] "output-XETG00089__0011236__11236X_5459B_HYP__20240223__174831"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 36415037
# [1] 30282008
# [1] 0.8315798
# [1] "output-XETG00089__0011236__11236X_8667C_HYP__20240223__174831"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 30585232
# [1] 23950376
# [1] 0.7830699


##################################################################################################
######## filter out low-quality transcripts and create csv files for within-cell and with-nucleus data
##################################################################################################
# library(fread)
samples_filderName = list.files(paste0(dir_data))
for(sample_tmp in samples_filderName[7]){
  print(sample_tmp)
  ## read transcripts.csv 
  # gunzip -k transcripts.csv.gz
  df_transcripts <- fread(paste0(dir_data,sample_tmp,
                                 # "/output-XETG00089__00011586__Region_3__20231019__161214/",
                                 "/transcripts.csv.gz"), header = T)
  ###### filter the transcripts
  df_transcripts_QC20 = df_transcripts[df_transcripts$qv >= 20,]
  
  ###### create a csv file with only filtered transcripts
  dir.create(paste0(dir_data_processed,"/cellLevel/",sample_tmp,"_cellLevel_qv20"))
  
  print("writing csv file")
  
  fwrite(df_transcripts_QC20, file= paste0(dir_data_processed,"/cellLevel/",sample_tmp,"_cellLevel_qv20/", "transcripts.csv"))
  
  ###### keep with-nuclei data
  dir.create(paste0(dir_data_processed,"/nucleusLevel/",sample_tmp,"_nucleusLevel_qv20"))

  df_transcripts_QC20_nuclear = df_transcripts_QC20[df_transcripts_QC20$overlaps_nucleus ==1  ,]
  fwrite(df_transcripts_QC20_nuclear, file= paste0(dir_data_processed,"/nucleusLevel/",sample_tmp,"_nucleusLevel_qv20/", "transcripts.csv"))
  
}

unique(df_transcripts$overlaps_nucleus)
# [1] 0 1
# 
head(df_transcripts)
# transcript_id    cell_id overlaps_nucleus feature_name x_location
# <i64>     <char>            <int>       <char>      <num>
#   1: 281814279127098 UNASSIGNED                0        CELF2   13.54553
# 2: 281814279127152 UNASSIGNED                0        BMP8B   27.48775
# 3: 281814279127183 UNASSIGNED                0        CELF2   34.71953
# 4: 281814279127237 UNASSIGNED                0        CELF2   48.32832
# 5: 281814279127532 UNASSIGNED                0        CELF2   90.61813
# 6: 281814279127589 UNASSIGNED                0        CELF2   95.82365
# y_location z_location       qv fov_name nucleus_distance
# <num>      <num>    <num>   <char>            <num>
#   1:   5943.480   13.05608 40.00000      AG1         324.9864
# 2:   5911.152   13.20928 16.23016      AG1         350.3790
# 3:   5997.533   13.01907 40.00000      AG1         266.9502
# 4:   5940.587   12.99224 40.00000      AG1         316.0282
# 5:   5768.442   13.03705 40.00000      AG1         442.1286
# 6:   5890.335   13.15292 40.00000      AG1         338.8733

dim((df_transcripts))
# [1] 30,519,510       10


##################################################################################################
##### copy and paste the files to a organized folder set, grouped by sample
##################################################################################################
# unchanged data
for(sample_tmp in samples_filderName){
  print(sample_tmp)
  
  files_nucleo = c(paste0(dir_data,sample_tmp, "/nucleus_boundaries.csv.gz"),
                      paste0(dir_data,sample_tmp, "/nucleus_boundaries.parquet") )
  files_cell = c(paste0(dir_data,sample_tmp, "/cell_boundaries.csv.gz"),
                      paste0(dir_data,sample_tmp, "/cell_boundaries.parquet"))
  
  ###### cell-level data
  file.copy(files_cell, paste0(dir_data_processed,"/cellLevel/",sample_tmp,"_cellLevel_qv20/"))

  ###### nucleus-level data
  file.copy(files_nucleo, paste0(dir_data_processed,"/nucleusLevel/",sample_tmp,"_nucleusLevel_qv20/"))
 
}








