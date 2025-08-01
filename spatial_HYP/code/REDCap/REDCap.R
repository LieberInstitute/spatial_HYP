setwd('/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/')

library("here")

REDCap <- read.csv(file.path(here::here("raw-data", "sample_info_Visium", "Visium_DATA_2023-05-17_1500.csv")), header = TRUE, stringsAsFactors = FALSE)
A1 <- subset(REDCap, select = c("slide", "species_a1", "sample_a1", "serial_a1","region_a1", "project_a1"))
B1 <- subset(REDCap, select = c("slide", "species_b1", "sample_b1", "serial_b1","region_b1", "project_b1"))
C1 <- subset(REDCap, select = c("slide", "species_c1", "sample_c1", "serial_c1","region_c1", "project_c1"))
D1 <- subset(REDCap, select = c("slide", "species_d1", "sample_d1", "serial_d1","region_d1", "project_d1"))
colnames(A1) <- colnames(B1) <- colnames(C1) <- colnames(D1) <- c("slide", "species", "sample", "serial","region", "project")
A1$array <- "A1"
B1$array <- "B1"
C1$array <- "C1"
D1$array <- "D1"

REDCap_table <- rbind(A1, B1, C1, D1)
REDCap_table <- REDCap_table[order(REDCap_table$slide), ]

REDCap_HYP <- REDCap_table[which(REDCap_table$project == "LIBD4195"), ]

Brain_nums <- unique(REDCap_HYP$sample)
write.table(Brain_nums, file = (here::here("raw-data", "sample_info_Visium", "ALLbrains.txt")), row.names = FALSE, col.names = FALSE)

Samples <- unique(paste0(REDCap_HYP$slide, "_", REDCap_HYP$array))
write.table(Samples, file = (here::here("raw-data", "sample_info_Visium", "ALLsamples.txt")), row.names = FALSE, col.names = FALSE)

save(REDCap_HYP, file = (here::here("code", "REDCap", "REDCap_HYP.rda")))

## https://jhu-genomics.slack.com/archives/CR9NYA0BF/p1650383126365919
## /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2022-04-12_SPag033122/
# dir.create(here::here("raw-data", "FASTQ", "2022-04-12_SPag033122"))
# 
# samples <- subset(REDCap_HPC, select = c("date", "slide", "array", "sample_number"))
# samples$date <- as.Date(samples$date)
# samples <- samples[samples$date > "2021-10-11", ]
# samples$samples <- paste0(samples$slide, "_", samples$array)
# write.table(samples$samples, file = (here::here("code", "01_spaceranger", "samples_2022-04-12_SPag033122.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# samples$rawFASTQ <- paste0("/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2022-04-12_SPag033122/", samples$sample_number, "*")
# samples$softlink <- paste("raw-data", "FASTQ", "2022-04-12_SPag033122", samples$samples, sep = "/")
# sapply(samples$softlink, dir.create)
# 
# temp <- as.factor(paste0("ln -s ", samples$rawFASTQ, " ", samples$softlink, "/"))
# write.table(temp, file = (here::here("raw-data", "FASTQ", "2022-04-12_SPag033122", "README.md")), row.names = FALSE, col.names = FALSE, quote = FALSE)
