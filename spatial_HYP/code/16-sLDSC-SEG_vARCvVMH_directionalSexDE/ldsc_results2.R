
library(data.table)
# read data 
setwd("/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/spatial_HYP/code/16-sLDSC-SEG_vARCvVMH_directionalSexDE")
dat <- fread("ldsc_results.txt",as.is=T,header=T,sep="\t")

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
dat$trait[trait=="GSCAN_AgeSmk"] <- "Age of smoking"
dat$trait[trait=="GSCAN_CigDay"] <- "Cigarettes per day"
dat$trait[trait=="GSCAN_DrnkWk"] <- "Drinks per week"
dat$trait[trait=="GSCAN_SmkCes"] <- "Smoking cessation"
dat$trait[trait=="GSCAN_SmkInit"] <- "Smoking initiation"
dat$trait[trait=="Schizophrenia_PGC3"] <- "Schizophrenia (PGC3)"
dat$trait[trait=="Bipolar Disorder2"] <- "Bipolar (PGC2)"
dat$trait[trait=="epilepsy"] <- "Epilepsy"
dat$trait[trait=="Type_2_Diabetes"] <- "Type 2 Diabetes"
dat$trait[trait=="Alzheimer Disease3"] <- "Alzheimer Disease v3"
dat$trait[trait=="mdd2019edinburgh"] <- "Depression"

tmpdat <- copy(as.data.table(dat))

# FDR correction
dat$p_zcore <- pnorm(abs(dat$Coefficient_z.score),lower.tail=F)*2
dat$FDR <- p.adjust(dat$p_zcore,method="fdr")


write.csv(dat,"ldsc_results_directional_SexDEsetsFDR_vVMH_vARC_only.csv")