# ml conda_R/4.3.x
# R

###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/"
dir_data_processed = "/dcs04/hansen/data/ywang/xenium/data_processed/"
dir_out = "/dcs04/hansen/data/ywang/xenium/out/"

###########
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
###########
ColorOut()
options("styler.addins_style_transformer" = "biocthis::bioc_style()")

## change write.table to output a TSV without quotes or rownames by default
default(write.table) <- list(row.names=F,col.names=T,sep='\t',quote=F)

theme_set(theme_bw()+theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.y = element_text(size =16), plot.title = element_text(size = 20,hjust=0.5), strip.text = element_text(size=18), legend.text = element_text(size=10), legend.title = element_text(size=11,hjust=0.5)))





####### MOLECULEEXPERIMENT to SpatialExperiment #######

## load all Xenium results from their parent directory; change to understandable sample names; import cell (not nucleus) boundaries
xens <- MoleculeExperiment::readXenium(paste0(dir_data,"/20231019__161052__101923_KMon101923/"),keepCols = "essential",addBoundaries = "cell")

### rename samples (which are the names of the molecules and boundaries parts of the molecule experiment)
# sample/section IDs by donor, spatial reps as letters; in order of slide #86 region 1-2-3, #97 reg 1-2-3 (how they will be loaded by default)
xennames <- c("br6197","br5993a","br5993b","br1735a","br1735b","br6588")
i<-1
for (i in c(1:length(xens@molecules$detected))){
  names(xens@molecules$detected)[i] <- xennames[i]
  names(xens@boundaries$cell)[i] <- xennames[i]
}
ls(names(xens@molecules$detected)[[1]])
# > names(xens@molecules$detected[[1]])
# [1] "ABCC9"                   "ACVR1C"                 
# [3] "ADAMTS12"                "ADAMTS16"               
# [5] "ADAMTS3"                 "ADARB1"                 
# [7] "ADCYAP1"                 "ADRA1A"                 
# [9] "ADRA1B"                  "AGRP"                   
# [11] "AIF1"                    "ALK"                    
# [13] "ALKAL2"                  "ANGPT1"                 
# [15] "ANK1"                    "ANKRD18A"               
# [17] "ANKRD34B"                "ANO3"                   
# [19] "ANXA1"                   "APOE"                   
# [21] "APP"                     "AQP4"                   
# [23] "AR"                      "ARHGAP24"               
# [25] "ATP2C2"                  "AVP_gene"               
# [27] "B4GALNT1"                "BCAN"                   
# [29] "BDNF"                    "BLANK_0001"             
# [31] "BLANK_0005"              "BLANK_0007"             
# [33] "BLANK_0009"              "BLANK_0012"             
# [35] "BLANK_0013"              "BLANK_0015"             
# [37] "BLANK_0020"              "BLANK_0021"             
# [39] "BLANK_0023"              "BLANK_0024"             
# [41] "BLANK_0027"              "BLANK_0028"             
# [43] "BLANK_0031"              "BLANK_0036"             
# [45] "BLANK_0039"              "BLANK_0041"             
# [47] "BLANK_0043"              "BLANK_0044"             
# [49] "BLANK_0045"              "BLANK_0054"             
# [51] "BLANK_0055"              "BLANK_0056"             
# [53] "BLANK_0058"              "BLANK_0059"             
# [55] "BLANK_0062"              "BLANK_0063"             
# [57] "BLANK_0066"              "BLANK_0067"             
# [59] "BLANK_0068"              "BLANK_0071"             
# [61] "BLANK_0075"              "BLANK_0078"             
# [63] "BLANK_0079"              "BLANK_0080"             
# [65] "BLANK_0086"              "BLANK_0093"             
# [67] "BLANK_0095"              "BLANK_0097"             
# [69] "BLANK_0098"              "BLANK_0100"             
# [71] "BLANK_0103"              "BLANK_0109"             
# [73] "BLANK_0116"              "BLANK_0117"             
# [75] "BLANK_0118"              "BLANK_0124"             
# [77] "BLANK_0127"              "BLANK_0128"             
# [79] "BLANK_0129"              "BLANK_0136"             
# [81] "BLANK_0139"              "BLANK_0140"             
# [83] "BLANK_0143"              "BLANK_0146"             
# [85] "BLANK_0150"              "BLANK_0153"             
# [87] "BLANK_0156"              "BLANK_0161"             
# [89] "BLANK_0164"              "BLANK_0167"             
# [91] "BLANK_0169"              "BLANK_0172"             
# [93] "BLANK_0173"              "BLANK_0179"             
# [95] "BLANK_0185"              "BLANK_0190"             
# [97] "BLANK_0192"              "BLANK_0194"             
# [99] "BLANK_0195"              "BLANK_0205"             
# [101] "BLANK_0208"              "BLANK_0209"             
# [103] "BLANK_0213"              "BLANK_0217"             
# [105] "BLANK_0224"              "BLANK_0225"             
# [107] "BLANK_0231"              "BLANK_0233"             
# [109] "BLANK_0235"              "BLANK_0236"             
# [111] "BLANK_0238"              "BLANK_0245"             
# [113] "BLANK_0247"              "BLANK_0251"             
# [115] "BLANK_0252"              "BLANK_0254"             
# [117] "BLANK_0260"              "BLANK_0272"             
# [119] "BLANK_0274"              "BLANK_0277"             
# [121] "BLANK_0281"              "BLANK_0290"             
# [123] "BLANK_0291"              "BLANK_0293"             
# [125] "BLANK_0295"              "BLANK_0314"             
# [127] "BLANK_0317"              "BLANK_0321"             
# [129] "BLANK_0345"              "BLANK_0346"             
# [131] "BLANK_0366"              "BLANK_0368"             
# [133] "BLANK_0370"              "BLANK_0376"             
# [135] "BLANK_0381"              "BLANK_0393"             
# [137] "BMP8B"                   "BRINP3"                 
# [139] "BTBD11"                  "C1QL3"                  
# [141] "C1orf162"                "C2orf80"                
# [143] "CALCRL"                  "CAMK2A"                 
# [145] "CAPG"                    "CAPN3"                  
# [147] "CAV1"                    "CCK"                    
# [149] "CCL4"                    "CCL5"                   
# [151] "CCNA1"                   "CCNB2"                  
# [153] "CD14"                    "CD163"                  
# [155] "CD2"                     "CD36"                   
# [157] "CD3G"                    "CD4"                    
# [159] "CD48"                    "CD52"                   
# [161] "CD68"                    "CD83"                   
# [163] "CD86"                    "CDH1"                   
# [165] "CDH12"                   "CDH13"                  
# [167] "CDH4"                    "CDH6"                   
# [169] "CDK1"                    "CELF2"                  
# [171] "CEMIP"                   "CEMIP2"                 
# [173] "CENPF"                   "CHGB"                   
# [175] "CHODL"                   "CLDN11"                 
# [177] "CNDP1"                   "CNTN2"                  
# [179] "CNTNAP3B"                "COL12A1"                
# [181] "COL25A1"                 "CORO1A"                 
# [183] "CRH"                     "CRHBP"                  
# [185] "CRHR2"                   "CRYM"                   
# [187] "CSPG4"                   "CTNNA3"                 
# [189] "CTSH"                    "CTSS"                   
# [191] "CUX2"                    "CX3CR1"                 
# [193] "CXCL14"                  "CXCR4"                  
# [195] "CYP19A1"                 "CYP26A1"                
# [197] "CYP26B1"                 "CYTIP"                  
# [199] "DCC"                     "DCN"                    
# [201] "DDN"                     "DDR2"                   
# [203] "DLG4"                    "DNER"                   
# [205] "DOC2A"                   "DYDC2"                  
# [207] "DeprecatedCodeword_0160" "DeprecatedCodeword_0171"
# [209] "DeprecatedCodeword_0196" "DeprecatedCodeword_0304"
# [211] "DeprecatedCodeword_0318" "DeprecatedCodeword_0344"
# [213] "DeprecatedCodeword_0373" "ECEL1"                  
# [215] "EFHD1"                   "EGFR"                   
# [217] "ELK1"                    "ELOVL2"                 
# [219] "ENOX2"                   "ERBB3"                  
# [221] "ERMN"                    "ESR1"                   
# [223] "EYA4"                    "FABP6"                  
# [225] "FASLG"                   "FBLN1"                  
# [227] "FCER1G"                  "FCGR1A"                 
# [229] "FCGR3A"                  "FEZF1"                  
# [231] "FGFR2"                   "FGFR3"                  
# [233] "FILIP1"                  "FLT1"                   
# [235] "FOS"                     "FRMPD2"                 
# [237] "FSTL4"                   "GABRA5"                 
# [239] "GABRE"                   "GABRQ"                  
# [241] "GAD1"                    "GAD2"                   
# [243] "GAL"                     "GAS2L3"                 
# [245] "GDF10"                   "GFOD1"                  
# [247] "GHRH"                    "GJA1"                   
# [249] "GLRA2"                   "GLRA3"                  
# [251] "GNLY"                    "GPER1"                  
# [253] "GPNMB"                   "GPR183"                 
# [255] "GPR34"                   "GSX1"                   
# [257] "GZMA"                    "HCRT"                   
# [259] "HES1"                    "HHATL"                  
# [261] "HILPDA"                  "HLA-DMB"                
# [263] "HLA-DQA1"                "HRH1"                   
# [265] "HRH3"                    "HS3ST2"                 
# [267] "HS3ST4"                  "HS6ST2"                 
# [269] "HTR2A"                   "HTR2C"                  
# [271] "ICAM5"                   "IDH1"                   
# [273] "IDH2"                    "IDO1"                   
# [275] "IFITM3"                  "IGFBP3"                 
# [277] "IGFBP4"                  "IGFBP5"                 
# [279] "IL18R1"                  "IL7R"                   
# [281] "ITGA8"                   "ITGAM"                  
# [283] "ITGAX"                   "ITGB2"                  
# [285] "KCNH5"                   "KIT"                    
# [287] "KLF2"                    "KLF4"                   
# [289] "KLK6"                    "KLRB1"                  
# [291] "KRT18"                   "KRT8"                   
# [293] "LAMA2"                   "LAMP5"                  
# [295] "LBHD2"                   "LEPR"                   
# [297] "LHX6"                    "LINC01615"              
# [299] "LINC02763"               "LOX"                    
# [301] "LRRK1"                   "LRRK2"                  
# [303] "LY86"                    "LYPD6"                  
# [305] "LYPD6B"                  "LYVE1"                  
# [307] "MAG"                     "MAL"                    
# [309] "MAN1A1"                  "MC3R"                   
# [311] "MC4R"                    "MCTP2"                  
# [313] "MEIS2"                   "MEPE"                   
# [315] "MGST1"                   "MKI67"                  
# [317] "MOBP"                    "MOG"                    
# [319] "MS4A6A"                  "MT2A"                   
# [321] "MYO16"                   "MYO5B"                  
# [323] "MYRF"                    "MYT1L"                  
# [325] "MYZAP"                   "NCSTN"                  
# [327] "NDST4"                   "NECAB2"                 
# [329] "NES"                     "NKG7"                   
# [331] "NNAT"                    "NOTCH1"                 
# [333] "NPFFR2"                  "NPNT"                   
# [335] "NPY1R"                   "NPY2R"                  
# [337] "NR2F2"                   "NR4A2"                  
# [339] "NR5A1"                   "NRGN"                   
# [341] "NRP1"                    "NTNG1"                  
# [343] "NTNG2"                   "NTRK3"                  
# [345] "NWD2"                    "NXPH2"                  
# [347] "NegControlCodeword_0500" "NegControlCodeword_0501"
# [349] "NegControlCodeword_0502" "NegControlCodeword_0503"
# [351] "NegControlCodeword_0504" "NegControlCodeword_0505"
# [353] "NegControlCodeword_0506" "NegControlCodeword_0507"
# [355] "NegControlCodeword_0508" "NegControlCodeword_0509"
# [357] "NegControlCodeword_0510" "NegControlCodeword_0511"
# [359] "NegControlCodeword_0512" "NegControlCodeword_0513"
# [361] "NegControlCodeword_0514" "NegControlCodeword_0515"
# [363] "NegControlCodeword_0516" "NegControlCodeword_0517"
# [365] "NegControlCodeword_0518" "NegControlCodeword_0519"
# [367] "NegControlCodeword_0520" "NegControlCodeword_0521"
# [369] "NegControlCodeword_0522" "NegControlCodeword_0523"
# [371] "NegControlCodeword_0524" "NegControlCodeword_0525"
# [373] "NegControlCodeword_0526" "NegControlCodeword_0527"
# [375] "NegControlCodeword_0528" "NegControlCodeword_0529"
# [377] "NegControlCodeword_0530" "NegControlCodeword_0531"
# [379] "NegControlCodeword_0532" "NegControlCodeword_0533"
# [381] "NegControlCodeword_0534" "NegControlCodeword_0535"
# [383] "NegControlCodeword_0536" "NegControlCodeword_0537"
# [385] "NegControlCodeword_0538" "NegControlCodeword_0539"
# [387] "NegControlCodeword_0540" "NegControlProbe_00002"  
# [389] "NegControlProbe_00003"   "NegControlProbe_00004"  
# [391] "NegControlProbe_00009"   "NegControlProbe_00012"  
# [393] "NegControlProbe_00013"   "NegControlProbe_00014"  
# [395] "NegControlProbe_00016"   "NegControlProbe_00017"  
# [397] "NegControlProbe_00019"   "NegControlProbe_00022"  
# [399] "NegControlProbe_00024"   "NegControlProbe_00025"  
# [401] "NegControlProbe_00031"   "NegControlProbe_00033"  
# [403] "NegControlProbe_00034"   "NegControlProbe_00035"  
# [405] "NegControlProbe_00039"   "NegControlProbe_00041"  
# [407] "NegControlProbe_00042"   "OLIG1"                  
# [409] "OLIG2"                   "OPALIN"                 
# [411] "OTOGL"                   "OXTR"                   
# [413] "OXT_gene"                "P2RY12"                 
# [415] "P2RY13"                  "PAX6"                   
# [417] "PCDH19"                  "PCDH7"                  
# [419] "PCLO"                    "PCNA"                   
# [421] "PCP4"                    "PCSK1"                  
# [423] "PCSK6"                   "PDGFD"                  
# [425] "PDGFRA"                  "PECAM1"                 
# [427] "PEG10"                   "PHLDB2"                 
# [429] "PLCE1"                   "PLCH1"                  
# [431] "PLD5"                    "POMC"                   
# [433] "POSTN"                   "POU6F2"                 
# [435] "PROKR1"                  "PROX1"                  
# [437] "PSEN1"                   "PSEN2"                  
# [439] "PSENEN"                  "PTCHD4"                 
# [441] "PTGDS"                   "PTPRC"                  
# [443] "PTPRN"                   "PTPRZ1"                 
# [445] "PVALB"                   "PWP2"                   
# [447] "RAB3B"                   "RASGRP1"                
# [449] "RBFOX1"                  "RDH10"                  
# [451] "RELN"                    "RFTN1"                  
# [453] "RGS10"                   "RGS16"                  
# [455] "RIT2"                    "RNASET2"                
# [457] "RNF144B"                 "RNF152"                 
# [459] "RORB"                    "ROS1"                   
# [461] "RRAGB"                   "RSPO2"                  
# [463] "RXFP1"                   "RXRG"                   
# [465] "RYR3"                    "S100A4"                 
# [467] "SAMD5"                   "SDK1"                   
# [469] "SEMA5A"                  "SERPINA3"               
# [471] "SFRP2"                   "SHANK3"                 
# [473] "SHISAL2B"                "SLC10A4"                
# [475] "SLC17A6"                 "SLC17A7"                
# [477] "SLC24A3"                 "SLC26A4"                
# [479] "SLC6A3"                  "SLIT3"                  
# [481] "SNCA"                    "SNCG"                   
# [483] "SNTB2"                   "SORCS1"                 
# [485] "SOX10"                   "SOX11"                  
# [487] "SOX2"                    "SOX4"                   
# [489] "SOX9"                    "SPHKAP"                 
# [491] "SPI1"                    "SPON1"                  
# [493] "SST"                     "SSTR1"                  
# [495] "SSTR2"                   "ST18"                   
# [497] "STAT3"                   "STK32B"                 
# [499] "STXBP2"                  "SULF1"                  
# [501] "SYNPR"                   "TAC1"                   
# [503] "TAC3"                    "TACR1"                  
# [505] "TENM1"                   "TESPA1"                 
# [507] "TGFB1"                   "TGFB2"                  
# [509] "TGFBI"                   "THBS1"                  
# [511] "THEMIS"                  "THSD4"                  
# [513] "THSD7B"                  "TLL2"                   
# [515] "TMEM132C"                "TMEM233"                
# [517] "TMIGD3"                  "TNRC6B"                 
# [519] "TOP2A"                   "TP53"                   
# [521] "TPH2"                    "TRAC"                   
# [523] "TREM2"                   "TRH"                    
# [525] "TRHDE"                   "TRIL"                   
# [527] "TRPC5"                   "TRPC6"                  
# [529] "TSHZ2"                   "TTYH1"                  
# [531] "UGT8"                    "UNC5B"                  
# [533] "VCAN"                    "VGF"                    
# [535] "VIP"                     "VSNL1"                  
# [537] "VWC2L"                   "WIF1"                   
# [539] "ZBBX"                    "ZCCHC12"                
# [541] "ZDHHC23"                

### get names of features we want to keep (i.e. real genes)
keepgenes <- unlist(lapply(xens@molecules$detected,FUN=names))
keepgenes <- unique(keepgenes)
keepgenes <- grep(keepgenes,pattern="NegControlProbe|NegControlCodeword|DeprecatedCodeword|BLANK",value = T,invert = T)

### test gene subsetting on two samples
### subsetting MoleculeExperiment samples here doesn't have any straightforward approach off the top of my head (no function for it, anyhow), but i create a test object of just the first two samples to test the gene filtering. 
### ends up being a useful **demonstration of how to create a new, identical moleculeexperimentobject with only samples of interest.**
detectest <- list(xens@molecules$detected[[1]],xens@molecules$detected[[1]])
boundtest <- list(xens@boundaries$cell[[1]],xens@boundaries$cell[[2]])

names(detectest) <- xennames[c(1:2)]
names(boundtest) <- xennames[c(1:2)]
detectestholder <- list(detectest)
names(detectestholder) <- "detected"
boundtestholder <- list(boundtest)
names(boundtestholder) <- "cell"

test <- MoleculeExperiment(molecules=detectestholder,boundaries = boundtestholder)

i <-1
for (i in c(1:2)){
  test@molecules$detected[[i]] %<>% .[which(names(.) %in% keepgenes)]
}
stopifnot(length(names(test@molecules$detected[[1]]))==length(keepgenes))
## beautiful
rm(test,detectest,detectestholder,boundtest,boundtestholder,i)
gc(full=T)
### now really subset the moleculeexperiment
i<-1
for (i in c(1:6)){
  xens@molecules$detected[[i]] %<>% .[which(names(.) %in% keepgenes)]
}
rm(i)
gc(full=T)

### create a spatialexperiment for cell-level data ### 
### with the molecularexperiment function **countMolecules**. specify cell boundaries as the assay we want to use. (This is infinitely more straightforward, and looks to be the most computationally efficient route I've identified, for building an SPE, including rather than trying to build a SPE manually/piecewise from Xenium output files).

# This is able to perform multicore as long as BiocParallel is loaded.
xens <- countMolecules(xens,moleculesAssay = "detected",boundariesAssay = "cell",nCores = 8)

### holy SHIT that was fast. Like, 3 minutes.
gc(full=T)

# before we go further, we need to change the names of the "x_location" and "y_location" variables in the SPE's coldata, or ggspavis will get confused by that being names in both coldata and spatialcoords.
colnames(colData(xens))[3] <- "x_position"
colnames(colData(xens))[4] <- "y_position"

# save the resulting obj in the background
# job({saveRDS(xens,paste0(dir_data_processed, "hypxenia.cells.as.SPE_unfiltered.RDS"))#####start  from here!!!
  # job::export("none")},import="auto")
saveRDS(xens, file=paste0(dir_data_processed, "xens.rds"))
