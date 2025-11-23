
library(data.table)
# read data 
setwd("/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/spatial_HYP/code/16-sLDSC-SEG_vARCvVMH_directionalSexDE")
dat <- fread("ldsc_results.txt",header=T,sep='\t')

traits <- c(
"ADHD",
"Drinks Per Week",
"Alzheimer Disease",
"Alzheimer Disease2",
"Alzheimer Disease3",
"Anorexia",
"AUD_NatMed2023",
"Autism",
"BMI",
"Bipolar Disorder2",
"Bipolar Disorder",
"Alcohol Use Disorder",
"epilepsyAll",
"epilepsyFocal",
"epilepsy",
"GSCAN_AgeSmk",
"GSCAN_CigDay",
"GSCAN_DrnkWk",
"GSCAN_SmkCes",
"GSCAN_SmkInit",
"Height",
"Insomnia",
"Intelligence",
"mdd2019edinburgh",
"Depression_ex23andMe",
"Depression",
"STROKE1",
"STROKE2",
"OUD_MVP1",
"OUD_MVP12",
"OUD_MVP2",
"Type_2_Diabetes",
"Parkinson Disease",
"PTSD",
"Schizophrenia",
"Schizophrenia_PGC3",
"Cigarettes Per Day",
"STROKEall",
"Education Years",
"Neuroticism"
)

dat <- dat[trait %in% traits]

# rename traits
dat[trait=="GSCAN_AgeSmk",trait:="Age of smoking"]
dat[trait=="GSCAN_CigDay",trait:="Cigarettes per day"]
dat[trait=="GSCAN_DrnkWk",trait:="Drinks per week"]
dat[trait=="GSCAN_SmkCes",trait:="Smoking cessation"]
dat[trait=="GSCAN_SmkInit",trait:="Smoking initiation"]
dat[trait=="Schizophrenia_PGC3",trait:="Schizophrenia (PGC3)"]
dat[trait=="Bipolar Disorder2",trait:="Bipolar (PGC2)"]
dat[trait=="epilepsy",trait:="Epilepsy"]
dat[trait=="Type_2_Diabetes",trait:="Type 2 Diabetes"]
dat[trait=="Alzheimer Disease3",trait:="Alzheimer Disease v3"]
dat[trait=="mdd2019edinburgh",trait:="Depression"]

tmpdat <- copy(as.data.table(dat))

# FDR correction
dat[,p_zcore:=pnorm(abs(`Coefficient_z-score`),lower.tail=F)*2]
dat[,FDR:=p.adjust(dat$p_zcore,method="fdr")]


write.csv(dat,"ldsc_results_directional_SexDEsetsFDR_vVMH_vARC_only.csv")
