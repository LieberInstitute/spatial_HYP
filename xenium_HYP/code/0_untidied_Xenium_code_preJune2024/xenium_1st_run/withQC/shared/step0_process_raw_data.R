# ml conda_R/4.3.x
# R
####################################
####################################
####################################
####################################
####################################
####################################
library(data.table)
library(Biostrings)
library(ggplot2)
library(gridExtra)
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

####################################
####################################
####################################
cutoff_QC =  20
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20231019__161052__101923_KMon101923/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed_QC/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out_QC/"

####################################
####################################
dir.create(paste0(dir_data_processed,"/cellLevel/"))
dir.create(paste0(dir_data_processed,"/nucleusLevel/"))
####################################
####################################
samples_filderName = list.files(paste0(dir_data))
for(sample_tmp in samples_filderName){
  print(sample_tmp)
  ## read transcripts.csv 
  # gunzip -k transcripts.csv.gz
  df_transcripts <- fread(paste0(dir_data,sample_tmp,
                                 # "/output-XETG00089__00011586__Region_3__20231019__161214/",
                                 "/transcripts.csv.gz"), header = T)
  ###### filter the transcripts
  df_transcripts_QC20 = df_transcripts[df_transcripts$qv >= cutoff_QC,]
  
  ###### create a csv file with only filtered transcripts
  dir.create(paste0(dir_data_processed,"/cellLevel/",sample_tmp,"_cellLevel_qv20"))
  
  print("writing csv file")
  
  fwrite(df_transcripts_QC20, file= paste0(dir_data_processed,"/cellLevel/",sample_tmp,"_cellLevel_qv20/", "transcripts.csv"))
  
  ###### keep with-nuclei data
  dir.create(paste0(dir_data_processed,"/nucleusLevel/",sample_tmp,"_nucleusLevel_qv20"))

  df_transcripts_QC20_nuclear = df_transcripts_QC20[df_transcripts_QC20$overlaps_nucleus ==1  ,]
  fwrite(df_transcripts_QC20_nuclear, file= paste0(dir_data_processed,"/nucleusLevel/",sample_tmp,"_nucleusLevel_qv20/", "transcripts.csv"))
  
}


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

##################################################################################################
### nucleus level processing
##################################################################################################
######## laod raw data
xens <- MoleculeExperiment::readXenium(paste0(dir_data_processed,"/nucleusLevel/"),
                                       keepCols = "essential",
                                       addBoundaries = "nucleus")  


######## 
for (i in c(1:samples_filderName)){
  names(xens@molecules$detected)[i] <- samples_filderName[i]
  names(xens@boundaries$nucleus)[i] <- samples_filderName[i]
}

######## this step takes ~30min for running 13 samples
xens <- countMolecules(xens,moleculesAssay = "detected",boundariesAssay = "nucleus")

######## 
return(xens)
