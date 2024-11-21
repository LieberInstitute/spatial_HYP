setwd("/Users/atharv.chandra/Desktop/HYP_1225")
#setwd("/Volumes/neural_plasticity/drd5_projects/Human/NAc/new images/6423")

#Load some more libraries if needed
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(pracma)
library(RColorBrewer)
library(writexl)
library(readxl)
library(stringr)


#copy the data
data_copy <- combined_data

#set conserved variables & paths
threshold <- 0.5
width_img <- (16/3)
height_img <- (18/3)
pdf_title_whole <- "Left & Right HYP Slices Whole Three Clusters"
pdf_title_polygon <- "Left & Right HYP Slices Polygon Three Clusters"
pdf_title_selected <- "Left & Right HYP Slices Selected Three Clusters"
pre_file_path_single_channels <- "/Users/atharv.chandra/Desktop/HYP_1225/Single Channel plots"
pre_file_path_selecteds <- "/Users/atharv.chandra/Desktop/HYP_1225/Selections"
pre_file_path_results <- "/Users/atharv.chandra/Desktop/HYP_1225"

#find all positive cells by channel
pos_520 <- combined_data[kmeans_result_selected_column_data_520$cluster %in% top_two_clusters_520, ]
pos_570 <- combined_data[kmeans_result_selected_column_data_570$cluster %in% top_two_clusters_570, ]
pos_690 <- combined_data[kmeans_result_selected_column_data_690$cluster %in% top_two_clusters_690, ]

#Store the object IDs of each as strings (for later phenotyping)
pos_520$row_number <- rownames(pos_520)
pos_570$row_number <- rownames(pos_570)
pos_690$row_number <- rownames(pos_690)

#Plotting

#Tell script where to save the file and under what name
file_path <- file.path(pre_file_path_single_channels, paste(pdf_title_whole, titles[1], ".pdf"))
#Make a pdf canvas for the script to plot on, with dimensions specified earlier
pdf(file = file_path, width = width_img, height = height_img)

#save 520 plot by channel
plot(data_copy$XMax, data_copy$YMax, 
     xlim = c(0, max(data_copy$XMax)), 
     ylim = c(0, max(data_copy$YMax)), 
     col = "lightgrey", pch = 20, cex = .2, title(paste(pdf_title_whole, titles[1])))
points(pos_520$XMax, pos_520$YMax, col = "green", pch = 20, cex = .2)
dev.off()

#save 570 plot by channel
file_path <- file.path(pre_file_path_single_channels, paste(pdf_title_whole, titles[2], ".pdf"))
pdf(file = file_path, width = width_img, height = height_img)
plot(data_copy$XMax, data_copy$YMax, 
     xlim = c(0, max(data_copy$XMax)), 
     ylim = c(0, max(data_copy$YMax)), 
     col = "lightgrey", pch = 20, cex = .2, title(paste(pdf_title_whole, titles[2])))
points(pos_570$XMax, pos_570$YMax, col = "blue", pch = 20, cex = .2)
dev.off()

#save 620 plot by channel
file_path <- file.path(pre_file_path_single_channels, paste(pdf_title_whole, titles[3], ".pdf"))
pdf(file = file_path, width = width_img, height = height_img)
plot(data_copy$XMax, 
     data_copy$YMax, 
     xlim = c(0, max(data_copy$XMax)), 
     ylim = c(0, max(data_copy$YMax)), 
     col = "lightgrey", pch = 20, cex = .2, title(paste(pdf_title_whole, titles[3])))
points(pos_690$XMax, pos_690$YMax, col = "purple", pch = 20, cex = .2)
dev.off()

#Import polygons

#pull the data from the csv
polygon1_points_read <- read_excel(path = file.path(pre_file_path_selecteds, paste(pdf_title_selected,"polys", ".xlsx")), sheet = "Polygon 1 Points")
#check what cells are in the polygon
inside_polygon1 <- inpolygon(data_copy$XMax, data_copy$YMax, polygon1_points_read$x, polygon1_points_read$y)
#get object IDs
ingroup_ids1 <- data_copy$Object.Id[inside_polygon1]
#save the object ids in the same output form as the ROI function
poly1 <- data.frame(ingroup_ids1 = rep(NA_integer_, length(ingroup_ids1)))
poly1$ingroup_ids1 <- ingroup_ids1

#same but with the other polygon
polygon2_points_read <- read_excel(path = file.path(pre_file_path_selecteds, paste("03-",pdf_title_selected,"polys", ".xlsx")), sheet = "Polygon 2 Points")
#shifting the coords like the data was shifted
inside_polygon2 <- inpolygon(data_copy$XMax, data_copy$YMax, polygon2_points_read$x, polygon2_points_read$y)
ingroup_ids2 <- data_copy$Object.Id[inside_polygon2]
poly2 <- data.frame(ingroup_ids2 = rep(NA_integer_, length(ingroup_ids2)))
poly2$ingroup_ids2 <- ingroup_ids2

#Draw the polygons on the combined data
plot(data_copy$XMax, 
     data_copy$YMax, 
     xlim = c(0, max(data_copy$XMax)), 
     ylim = c(0, max(data_copy$YMax)), 
     col = "lightgrey", pch = 20, cex = .2, title(pdf_title_polygon))
polygon(polygon1_points_read$x, polygon1_points_read$y, 
        col = rgb(0, 0, 1), border = rgb(0, 0, 1))
polygon(polygon2_points_read$x, polygon2_points_read$y, 
        col = rgb(0, 0, 1), border = rgb(0, 0, 1))

# #Polygon Selector, uncomment if you don't have polys saved before
# create_polygons <- function(data) {
# 
#   #Make vectors to store x and y coordinates
#   polygon_x <- c()
#   polygon_y <- c()
#   ingroupids <- c()
#   #storing points Function
#   add_point <- function(x, y) {
#     polygon_x <<- c(polygon_x, x)
#     polygon_y <<- c(polygon_y, y)
#     points(x, y, col = "blue", pch = 20)
#   }
#   
#   plot(combined_data$XMax, combined_data$YMax, xlim = c(0, max(combined_data$XMax)), ylim = c(0, max(combined_data$YMax)), col = "lightgrey", pch = 20, cex = .2, title(paste("combined data", titles[2])))
#   #Click to add points to polygon
#   while (TRUE) {
#     #Get the coordinates of the mouse click
#     click <- locator(type = "p", n = 1)
# 
#     #If the user presses Escape, break the loop
#     if (is.null(click)) break
# 
#     #Add the clicked point to the polygon
#     add_point(click$x, click$y)
#   }
# 
#   #Draw the polygon
#   polygon(polygon_x, polygon_y, col = rgb(0, 0, 1), border = rgb(0, 0, 1))
# 
#   #Check if points are inside the polygon
#   inside_polygon <- inpolygon(data$XMax, data$YMax, polygon_x, polygon_y)
# 
#   #Plot random points, color them differently based on whether they are inside or outside the polygon
#   points(data$XMax[inside_polygon], data$YMax[inside_polygon], col = "green", pch = 20, cex = .2)
#   points(data$XMax[!inside_polygon], data$YMax[!inside_polygon], col = "red", pch = 20, cex = .2)
#   ingroupids <- data$Object.Id[inside_polygon]
# 
#   #get the output
#   result <- list(
#     polygon_points = data.frame(x = polygon_x, y = polygon_y),
#     ingroup_ids = ingroupids
#   )
# 
#   return(result)
# }
# #make a polygon
# poly1 <- create_polygons(data_copy)

# #make another polygon if needed
# poly2 <- create_polygons(data_copy)

# #save the polygon selections 
# write_xlsx(list("Polygon 1 Points" = poly1$polygon_points, "Polygon 2 Points" = poly2$polygon_points), path = file.path(pre_file_path_selecteds, paste(pdf_title_selected,"polys", ".xlsx")))

#combine selected in polygon cells
poly3 <- c(poly1$ingroup_ids1, poly2$ingroup_ids2)

#Add phenotype markers to data by checking if 
data_copy = data_copy %>%
  mutate(phenotype = case_when(
    (rownames(data_copy) %in% pos_520$row_number) & (rownames(data_copy) %in% pos_570$row_number) & (rownames(data_copy) %in% pos_690$row_number) ~ "520+570+690+",
    !(rownames(data_copy) %in% pos_520$row_number) & (rownames(data_copy) %in% pos_570$row_number) & (rownames(data_copy) %in% pos_690$row_number) ~ "520-570+690+",
    !(rownames(data_copy) %in% pos_520$row_number) & !(rownames(data_copy) %in% pos_570$row_number) & (rownames(data_copy) %in% pos_690$row_number) ~ "520-570-690+",
    (rownames(data_copy) %in% pos_520$row_number) & !(rownames(data_copy) %in% pos_570$row_number) & (rownames(data_copy) %in% pos_690$row_number)  ~ "520+570-690+",
    (rownames(data_copy) %in% pos_520$row_number) & (rownames(data_copy) %in% pos_570$row_number) & !(rownames(data_copy) %in% pos_690$row_number) ~ "520+570+690-",
    !(rownames(data_copy) %in% pos_520$row_number) & (rownames(data_copy) %in% pos_570$row_number) & !(rownames(data_copy) %in% pos_690$row_number) ~ "520-570+690-",
    !(rownames(data_copy) %in% pos_520$row_number) & !(rownames(data_copy) %in% pos_570$row_number) & !(rownames(data_copy) %in% pos_690$row_number) ~ "520-570-690-",
    (rownames(data_copy) %in% pos_520$row_number) & !(rownames(data_copy) %in% pos_570$row_number) & !(rownames(data_copy) %in% pos_690$row_number) ~ "520+570-690-"
  )) 

# change your data so it's just the stuff in the polygon
data_selected <- data_copy[poly3,]

#save a pdf of the selected data
file_path <- file.path(pre_file_path_selecteds, paste(pdf_title_polygon, ".pdf"))
pdf(file = file_path, width = width_img, height = height_img)
plot(data_copy$XMax, 
     data_copy$YMax, 
     xlim = c(0, max(data_copy$XMax)), 
     ylim = c(0, max(data_copy$YMax)), 
     col = "red", pch = 20, cex = .2, title(pdf_title_polygon))
points(data_selected$XMax, 
       data_selected$YMax, 
       col = "green", pch = 20, cex = .2)
dev.off()



#find all positive cells by channel of selected data by phenotype

#520
pos_520_selected <- data_selected %>%
  filter(str_detect(phenotype, "520\\+"))
#cell counts
length(pos_520_selected$Object.Id)

#570
pos_570_selected <- data_selected %>%
  filter(str_detect(phenotype, "570\\+"))
#cell counts
length(pos_570_selected$Object.Id)

#find all positive cells by channel of selected data by phenotype
pos_690_selected <- data_selected %>%
  filter(str_detect(phenotype, "690\\+"))
#cell counts
length(pos_690_selected$Object.Id)


# #save an image of selected area to your computer or server
file_path <- file.path(pre_file_path_selecteds, paste(pdf_title_selected, ".pdf"))
pdf(file = file_path, width = width_img, height = height_img)
plot(data_selected$XMax, 
     data_selected$YMax, 
     xlim = c(0, max(data_copy$XMax)), 
     ylim = c(0, max(data_copy$YMax)), 
     col = "lightgrey", pch = 20, cex = .2, title(pdf_title_selected))
points(pos_520_selected$XMax, pos_520_selected$YMax, col = "green", pch = 20, cex = .2)
points(pos_570_selected$XMax, pos_570_selected$YMax, col = "yellow", pch = 20, cex = .2)
points(pos_690_selected$XMax, pos_690_selected$YMax, col = "purple", pch = 20, cex = .2)
dev.off()

#save all data from all cells in selection
write_xlsx(data_selected, path = file.path(pre_file_path_results, paste("03-",pdf_title_selected, titles[1], titles[2], titles[3], ".xlsx")))

