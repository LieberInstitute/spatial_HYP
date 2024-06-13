# ml conda_R/4.3.x
###########
dir_data = "/dcs04/lieber/marmaypag/xeniumHYP_LIBD4195/xenium_hyp/raw-data/xenium/20240223__174743__022324_KMon120823/"
# dir.create("/dcs04/hansen/data/ywang/xenium_newData")

dir_data_processed = "/dcs04/hansen/data/ywang/xenium_newData/data_processed/"
dir_out = "/dcs04/hansen/data/ywang/xenium_newData/out/"

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


### Example seurat handling of the xenium data (up to the point of clustering, without any gene filters), is at the end of this Rmd. ###



####### MOLECULEEXPERIMENT to SpatialExperiment #######
# nucleus
## load all Xenium results from their parent directory; change to understandable sample names; import cell (not nucleus) boundaries
xens <- MoleculeExperiment::readXenium(paste0(dir_data,"/"),keepCols = "essential",addBoundaries = "cell")
# addBoundaries = "nucleus"


### rename samples (which are the names of the molecules and boundaries parts of the molecule experiment)
# sample/section IDs by donor, spatial reps as letters; in order of slide #86 region 1-2-3, #97 reg 1-2-3 (how they will be loaded by default)
xennames <- c("br1225a","br1225b","br8741c",
              "br8741d","br5459a","br5459b",
              "br8667c")
sex = c("M", "M","F", 
        "F", "M", "M",
        "F")
# 1225: Male
# 8667: Female
# 8741: Female
# 5459: Male

i<-1
for (i in c(1:length(xens@molecules$detected))){
  names(xens@molecules$detected)[i] <- xennames[i]
  names(xens@boundaries$cell)[i] <- xennames[i]
}
# > ls(xens@molecules$detected[[1]])
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
# [29] "BDNF"                    "BMP8B"                  
# [31] "BRINP3"                  "BTBD11"                 
# [33] "C1orf162"                "C1QL3"                  
# [35] "C2orf80"                 "CALCRL"                 
# [37] "CAMK2A"                  "CAPG"                   
# [39] "CAPN3"                   "CAV1"                   
# [41] "CCK"                     "CCL4"                   
# [43] "CCL5"                    "CCNA1"                  
# [45] "CCNB2"                   "CD14"                   
# [47] "CD163"                   "CD2"                    
# [49] "CD36"                    "CD3G"                   
# [51] "CD4"                     "CD48"                   
# [53] "CD52"                    "CD68"                   
# [55] "CD83"                    "CD86"                   
# [57] "CDH1"                    "CDH12"                  
# [59] "CDH13"                   "CDH4"                   
# [61] "CDH6"                    "CDK1"                   
# [63] "CELF2"                   "CEMIP"                  
# [65] "CEMIP2"                  "CENPF"                  
# [67] "CHGB"                    "CHODL"                  
# [69] "CLDN11"                  "CNDP1"                  
# [71] "CNTN2"                   "CNTNAP3B"               
# [73] "COL12A1"                 "COL25A1"                
# [75] "CORO1A"                  "CRH"                    
# [77] "CRHBP"                   "CRHR2"                  
# [79] "CRYM"                    "CSPG4"                  
# [81] "CTNNA3"                  "CTSH"                   
# [83] "CTSS"                    "CUX2"                   
# [85] "CX3CR1"                  "CXCL14"                 
# [87] "CXCR4"                   "CYP19A1"                
# [89] "CYP26A1"                 "CYP26B1"                
# [91] "CYTIP"                   "DCC"                    
# [93] "DCN"                     "DDN"                    
# [95] "DDR2"                    "DeprecatedCodeword_0160"
# [97] "DeprecatedCodeword_0171" "DeprecatedCodeword_0196"
# [99] "DeprecatedCodeword_0304" "DeprecatedCodeword_0318"
# [101] "DeprecatedCodeword_0344" "DeprecatedCodeword_0373"
# [103] "DLG4"                    "DNER"                   
# [105] "DOC2A"                   "DYDC2"                  
# [107] "ECEL1"                   "EFHD1"                  
# [109] "EGFR"                    "ELK1"                   
# [111] "ELOVL2"                  "ENOX2"                  
# [113] "ERBB3"                   "ERMN"                   
# [115] "ESR1"                    "EYA4"                   
# [117] "FABP6"                   "FASLG"                  
# [119] "FBLN1"                   "FCER1G"                 
# [121] "FCGR1A"                  "FCGR3A"                 
# [123] "FEZF1"                   "FGFR2"                  
# [125] "FGFR3"                   "FILIP1"                 
# [127] "FLT1"                    "FOS"                    
# [129] "FRMPD2"                  "FSTL4"                  
# [131] "GABRA5"                  "GABRE"                  
# [133] "GABRQ"                   "GAD1"                   
# [135] "GAD2"                    "GAL"                    
# [137] "GAS2L3"                  "GDF10"                  
# [139] "GFOD1"                   "GHRH"                   
# [141] "GJA1"                    "GLRA2"                  
# [143] "GLRA3"                   "GNLY"                   
# [145] "GPER1"                   "GPNMB"                  
# [147] "GPR183"                  "GPR34"                  
# [149] "GSX1"                    "GZMA"                   
# [151] "HCRT"                    "HES1"                   
# [153] "HHATL"                   "HILPDA"                 
# [155] "HLA-DMB"                 "HLA-DQA1"               
# [157] "HRH1"                    "HRH3"                   
# [159] "HS3ST2"                  "HS3ST4"                 
# [161] "HS6ST2"                  "HTR2A"                  
# [163] "HTR2C"                   "ICAM5"                  
# [165] "IDH1"                    "IDH2"                   
# [167] "IDO1"                    "IFITM3"                 
# [169] "IGFBP3"                  "IGFBP4"                 
# [171] "IGFBP5"                  "IL18R1"                 
# [173] "IL7R"                    "ITGA8"                  
# [175] "ITGAM"                   "ITGAX"                  
# [177] "ITGB2"                   "KCNH5"                  
# [179] "KIT"                     "KLF2"                   
# [181] "KLF4"                    "KLK6"                   
# [183] "KLRB1"                   "KRT18"                  
# [185] "KRT8"                    "LAMA2"                  
# [187] "LAMP5"                   "LBHD2"                  
# [189] "LEPR"                    "LHX6"                   
# [191] "LINC01615"               "LINC02763"              
# [193] "LOX"                     "LRRK1"                  
# [195] "LRRK2"                   "LY86"                   
# [197] "LYPD6"                   "LYPD6B"                 
# [199] "LYVE1"                   "MAG"                    
# [201] "MAL"                     "MAN1A1"                 
# [203] "MC3R"                    "MC4R"                   
# [205] "MCTP2"                   "MEIS2"                  
# [207] "MEPE"                    "MGST1"                  
# [209] "MKI67"                   "MOBP"                   
# [211] "MOG"                     "MS4A6A"                 
# [213] "MT2A"                    "MYO16"                  
# [215] "MYO5B"                   "MYRF"                   
# [217] "MYT1L"                   "MYZAP"                  
# [219] "NCSTN"                   "NDST4"                  
# [221] "NECAB2"                  "NegControlCodeword_0500"
# [223] "NegControlCodeword_0501" "NegControlCodeword_0502"
# [225] "NegControlCodeword_0503" "NegControlCodeword_0504"
# [227] "NegControlCodeword_0505" "NegControlCodeword_0506"
# [229] "NegControlCodeword_0507" "NegControlCodeword_0508"
# [231] "NegControlCodeword_0509" "NegControlCodeword_0510"
# [233] "NegControlCodeword_0511" "NegControlCodeword_0512"
# [235] "NegControlCodeword_0513" "NegControlCodeword_0514"
# [237] "NegControlCodeword_0515" "NegControlCodeword_0516"
# [239] "NegControlCodeword_0517" "NegControlCodeword_0518"
# [241] "NegControlCodeword_0519" "NegControlCodeword_0520"
# [243] "NegControlCodeword_0521" "NegControlCodeword_0522"
# [245] "NegControlCodeword_0523" "NegControlCodeword_0524"
# [247] "NegControlCodeword_0525" "NegControlCodeword_0526"
# [249] "NegControlCodeword_0527" "NegControlCodeword_0528"
# [251] "NegControlCodeword_0529" "NegControlCodeword_0530"
# [253] "NegControlCodeword_0531" "NegControlCodeword_0532"
# [255] "NegControlCodeword_0533" "NegControlCodeword_0534"
# [257] "NegControlCodeword_0535" "NegControlCodeword_0536"
# [259] "NegControlCodeword_0537" "NegControlCodeword_0538"
# [261] "NegControlCodeword_0539" "NegControlCodeword_0540"
# [263] "NegControlProbe_00002"   "NegControlProbe_00003"  
# [265] "NegControlProbe_00004"   "NegControlProbe_00009"  
# [267] "NegControlProbe_00012"   "NegControlProbe_00013"  
# [269] "NegControlProbe_00014"   "NegControlProbe_00016"  
# [271] "NegControlProbe_00017"   "NegControlProbe_00019"  
# [273] "NegControlProbe_00022"   "NegControlProbe_00024"  
# [275] "NegControlProbe_00025"   "NegControlProbe_00031"  
# [277] "NegControlProbe_00033"   "NegControlProbe_00034"  
# [279] "NegControlProbe_00035"   "NegControlProbe_00039"  
# [281] "NegControlProbe_00041"   "NegControlProbe_00042"  
# [283] "NES"                     "NKG7"                   
# [285] "NNAT"                    "NOTCH1"                 
# [287] "NPFFR2"                  "NPNT"                   
# [289] "NPY1R"                   "NPY2R"                  
# [291] "NR2F2"                   "NR4A2"                  
# [293] "NR5A1"                   "NRGN"                   
# [295] "NRP1"                    "NTNG1"                  
# [297] "NTNG2"                   "NTRK3"                  
# [299] "NWD2"                    "NXPH2"                  
# [301] "OLIG1"                   "OLIG2"                  
# [303] "OPALIN"                  "OTOGL"                  
# [305] "OXT_gene"                "OXTR"                   
# [307] "P2RY12"                  "P2RY13"                 
# [309] "PAX6"                    "PCDH19"                 
# [311] "PCDH7"                   "PCLO"                   
# [313] "PCNA"                    "PCP4"                   
# [315] "PCSK1"                   "PCSK6"                  
# [317] "PDGFD"                   "PDGFRA"                 
# [319] "PECAM1"                  "PEG10"                  
# [321] "PHLDB2"                  "PLCE1"                  
# [323] "PLCH1"                   "PLD5"                   
# [325] "POMC"                    "POSTN"                  
# [327] "POU6F2"                  "PROKR1"                 
# [329] "PROX1"                   "PSEN1"                  
# [331] "PSEN2"                   "PSENEN"                 
# [333] "PTCHD4"                  "PTGDS"                  
# [335] "PTPRC"                   "PTPRN"                  
# [337] "PTPRZ1"                  "PVALB"                  
# [339] "PWP2"                    "RAB3B"                  
# [341] "RASGRP1"                 "RBFOX1"                 
# [343] "RDH10"                   "RELN"                   
# [345] "RFTN1"                   "RGS10"                  
# [347] "RGS16"                   "RIT2"                   
# [349] "RNASET2"                 "RNF144B"                
# [351] "RNF152"                  "RORB"                   
# [353] "ROS1"                    "RRAGB"                  
# [355] "RSPO2"                   "RXFP1"                  
# [357] "RXRG"                    "RYR3"                   
# [359] "S100A4"                  "SAMD5"                  
# [361] "SDK1"                    "SEMA5A"                 
# [363] "SERPINA3"                "SFRP2"                  
# [365] "SHANK3"                  "SHISAL2B"               
# [367] "SLC10A4"                 "SLC17A6"                
# [369] "SLC17A7"                 "SLC24A3"                
# [371] "SLC26A4"                 "SLC6A3"                 
# [373] "SLIT3"                   "SNCA"                   
# [375] "SNCG"                    "SNTB2"                  
# [377] "SORCS1"                  "SOX10"                  
# [379] "SOX11"                   "SOX2"                   
# [381] "SOX4"                    "SOX9"                   
# [383] "SPHKAP"                  "SPI1"                   
# [385] "SPON1"                   "SST"                    
# [387] "SSTR1"                   "SSTR2"                  
# [389] "ST18"                    "STAT3"                  
# [391] "STK32B"                  "STXBP2"                 
# [393] "SULF1"                   "SYNPR"                  
# [395] "TAC1"                    "TAC3"                   
# [397] "TACR1"                   "TENM1"                  
# [399] "TESPA1"                  "TGFB1"                  
# [401] "TGFB2"                   "TGFBI"                  
# [403] "THBS1"                   "THEMIS"                 
# [405] "THSD4"                   "THSD7B"                 
# [407] "TLL2"                    "TMEM132C"               
# [409] "TMEM233"                 "TMIGD3"                 
# [411] "TNRC6B"                  "TOP2A"                  
# [413] "TP53"                    "TPH2"                   
# [415] "TRAC"                    "TREM2"                  
# [417] "TRH"                     "TRHDE"                  
# [419] "TRIL"                    "TRPC5"                  
# [421] "TRPC6"                   "TSHZ2"                  
# [423] "TTYH1"                   "UGT8"                   
# [425] "UnassignedCodeword_0001" "UnassignedCodeword_0005"
# [427] "UnassignedCodeword_0007" "UnassignedCodeword_0009"
# [429] "UnassignedCodeword_0012" "UnassignedCodeword_0013"
# [431] "UnassignedCodeword_0015" "UnassignedCodeword_0020"
# [433] "UnassignedCodeword_0021" "UnassignedCodeword_0023"
# [435] "UnassignedCodeword_0024" "UnassignedCodeword_0027"
# [437] "UnassignedCodeword_0028" "UnassignedCodeword_0031"
# [439] "UnassignedCodeword_0036" "UnassignedCodeword_0039"
# [441] "UnassignedCodeword_0041" "UnassignedCodeword_0043"
# [443] "UnassignedCodeword_0044" "UnassignedCodeword_0045"
# [445] "UnassignedCodeword_0054" "UnassignedCodeword_0055"
# [447] "UnassignedCodeword_0056" "UnassignedCodeword_0058"
# [449] "UnassignedCodeword_0059" "UnassignedCodeword_0062"
# [451] "UnassignedCodeword_0063" "UnassignedCodeword_0066"
# [453] "UnassignedCodeword_0067" "UnassignedCodeword_0068"
# [455] "UnassignedCodeword_0071" "UnassignedCodeword_0075"
# [457] "UnassignedCodeword_0078" "UnassignedCodeword_0079"
# [459] "UnassignedCodeword_0080" "UnassignedCodeword_0086"
# [461] "UnassignedCodeword_0093" "UnassignedCodeword_0095"
# [463] "UnassignedCodeword_0097" "UnassignedCodeword_0098"
# [465] "UnassignedCodeword_0100" "UnassignedCodeword_0103"
# [467] "UnassignedCodeword_0109" "UnassignedCodeword_0116"
# [469] "UnassignedCodeword_0117" "UnassignedCodeword_0118"
# [471] "UnassignedCodeword_0124" "UnassignedCodeword_0127"
# [473] "UnassignedCodeword_0128" "UnassignedCodeword_0129"
# [475] "UnassignedCodeword_0136" "UnassignedCodeword_0139"
# [477] "UnassignedCodeword_0140" "UnassignedCodeword_0143"
# [479] "UnassignedCodeword_0146" "UnassignedCodeword_0150"
# [481] "UnassignedCodeword_0153" "UnassignedCodeword_0156"
# [483] "UnassignedCodeword_0161" "UnassignedCodeword_0164"
# [485] "UnassignedCodeword_0167" "UnassignedCodeword_0169"
# [487] "UnassignedCodeword_0172" "UnassignedCodeword_0173"
# [489] "UnassignedCodeword_0179" "UnassignedCodeword_0185"
# [491] "UnassignedCodeword_0190" "UnassignedCodeword_0192"
# [493] "UnassignedCodeword_0194" "UnassignedCodeword_0195"
# [495] "UnassignedCodeword_0205" "UnassignedCodeword_0208"
# [497] "UnassignedCodeword_0209" "UnassignedCodeword_0213"
# [499] "UnassignedCodeword_0217" "UnassignedCodeword_0224"
# [501] "UnassignedCodeword_0225" "UnassignedCodeword_0231"
# [503] "UnassignedCodeword_0233" "UnassignedCodeword_0235"
# [505] "UnassignedCodeword_0236" "UnassignedCodeword_0238"
# [507] "UnassignedCodeword_0245" "UnassignedCodeword_0247"
# [509] "UnassignedCodeword_0251" "UnassignedCodeword_0252"
# [511] "UnassignedCodeword_0254" "UnassignedCodeword_0260"
# [513] "UnassignedCodeword_0272" "UnassignedCodeword_0274"
# [515] "UnassignedCodeword_0277" "UnassignedCodeword_0281"
# [517] "UnassignedCodeword_0290" "UnassignedCodeword_0291"
# [519] "UnassignedCodeword_0293" "UnassignedCodeword_0295"
# [521] "UnassignedCodeword_0314" "UnassignedCodeword_0317"
# [523] "UnassignedCodeword_0321" "UnassignedCodeword_0345"
# [525] "UnassignedCodeword_0346" "UnassignedCodeword_0366"
# [527] "UnassignedCodeword_0368" "UnassignedCodeword_0370"
# [529] "UnassignedCodeword_0376" "UnassignedCodeword_0381"
# [531] "UnassignedCodeword_0393" "UNC5B"                  
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
