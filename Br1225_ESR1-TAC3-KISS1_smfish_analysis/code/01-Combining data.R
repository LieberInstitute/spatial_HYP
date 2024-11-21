#Set working directory
setwd("/Users/atharv.chandra/Desktop/HYP_1225")

#import some libraries
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(pracma)
library(RColorBrewer)
library(writexl)
library(readxl)

#read the data 
data_left_name <- "Br1225_Left_HYP.2.csv"
data_right_name <- "Br1225_Right_HYP.2.csv"

data_left <- read.csv(data_left_name)
data_right <- read.csv(data_right_name)

#combine the data by shifting one slice to the right
data_right$XMax <- max(data_left$XMax)+data_right$XMax+500
data_right$XMin <- max(data_left$XMin)+data_right$XMin+500
data_right$YMax <- max(data_left$YMax)-max(data_right$YMax)+data_right$YMax
data_right$YMin <- max(data_left$YMin)-max(data_right$YMin)+data_right$YMin


data_right$Object.Id <- max(data_left$Object.Id)+data_right$Object.Id

#Then combine it into one dataframe
combined_data <- rbind(data_left, data_right)

#Flip the Y axis
combined_data$YMax <- max(combined_data$YMax)-combined_data$YMax
combined_data$YMin <- max(combined_data$YMin)-combined_data$YMin

#Plot to confirm
plot(combined_data$XMax, combined_data$YMax, xlim = c(0, max(combined_data$XMax)), ylim = c(0, max(combined_data$YMax)), col = "lightgrey", pch = 20, cex = .2, title("combined data"))


#Get titles and names 
columns <- c("X25xSil.Opal.520.Copies", "X20x.Opal.570.Copies", "X20x.Opal.690.Copies")

#gene titles
titles <- list(
  HYP = c("KISS1", "ESR1", "TAC3"),
  dlPFC = c("SLC17A7", "OPRM1", "GAD1", "DRD5"),
  Habenula = c("GPR151", "OPRM1", "CHAT", "DRD5"),
  NAc = c("PPP1R1B", "OPRM1", "CHAT", "DRD5"),
  dACC = c("PCP4", "DRD1", "VAT1L", "DRD5")
)

#relevant channel colors
colors <- c("green", "yellow", "purple")

#decide which title set you want to use
titles <- titles$HYP

