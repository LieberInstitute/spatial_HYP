#more libraries
library(factoextra)
library(cluster)

#select only columns of interest by channel
selected_520_data <- combined_data[, c("X25xSil.Opal.520.Copies", "X25xSil.Opal.520.Area..µm..")]
selected_690_data <- combined_data[, c("X20x.Opal.690.Copies", "X20x.Opal.690.Area..µm..")]
#cell intensity to confirm for Maddy
selected_570_data <- combined_data[, c("X20x.Opal.570.Copies", "X20x.Opal.570.Area..µm..")]

#Perform K means clustering by channel
kmeans_result_selected_column_data_570 <- kmeans(selected_570_data, centers = 3, nstart = 25)
kmeans_result_selected_column_data_520 <- kmeans(selected_520_data, centers = 3, nstart = 25)
kmeans_result_selected_column_data_690 <- kmeans(selected_690_data, centers = 3, nstart = 25)
#Display clusters
fviz_cluster(kmeans_result_selected_column_data_570, data = selected_570_data, geom = "point")
fviz_cluster(kmeans_result_selected_column_data_520, data = selected_520_data, geom = "point")
fviz_cluster(kmeans_result_selected_column_data_690, data = selected_690_data, geom = "point")

#Determine the cluster indices for the cell groups farthest from origin (0,0)

#690
distances_from_origin_690 <- 
  apply(kmeans_result_selected_column_data_690$centers, 1, function(center) {
  sqrt(center[1]^2 + center[2]^2)  # Euclidean distance formula: sqrt(X^2 + Y^2)
})
sorted_indices_690 <- order(distances_from_origin_690, decreasing = TRUE)
top_two_clusters_690 <- sorted_indices_690[1:2]

#520
distances_from_origin_520 <- 
  apply(kmeans_result_selected_column_data_520$centers, 1, function(center) {
    sqrt(center[1]^2 + center[2]^2)  # Euclidean distance formula: sqrt(X^2 + Y^2)
  })
sorted_indices_520 <- order(distances_from_origin_520, decreasing = TRUE)
top_two_clusters_520 <- sorted_indices_520[1:2]

#570
distances_from_origin_570 <- 
  apply(kmeans_result_selected_column_data_570$centers, 1, function(center) {
    sqrt(center[1]^2 + center[2]^2)  # Euclidean distance formula: sqrt(X^2 + Y^2)
  })
sorted_indices_570 <- order(distances_from_origin_570, decreasing = TRUE)
top_two_clusters_570 <- sorted_indices_570[1:2]


#Plot by channel & return number of cells by cluster
#PAY ATTENTION TO CLUSTER NUMBER & CLUSTER PAIRS! 
#CLUSTER 1 DOES NOT ALWAYS MEAN THE BOTTOM CLUSTER
#USE THE PLOTS BELOW TO FIND CORRECT CLUSTER

#Spatially Plot by Channel



# Sort the distances in descending order and get the indices of the two largest distances


# 2. Plot the spatial distribution of points from 'combined_data' for the two clusters
plot(combined_data$XMax, 
     combined_data$YMax, 
     xlim = c(0, max(combined_data$XMax)), 
     ylim = c(0, max(combined_data$YMax)), 
     col = "lightgrey", pch = 20, cex = .2, 
     main = "Spatial Plot for Top Two Clusters: 520 nm Channel")


# Loop through the top two clusters and spatially plot 
#   cells and print cell counts by cluster

#690
plot(combined_data$XMax, 
     combined_data$YMax, 
     xlim = c(0, max(combined_data$XMax)), 
     ylim = c(0, max(combined_data$YMax)), 
     col = "lightgrey", pch = 20, cex = .2, title(paste(titles[3])))
for (cluster_id in top_two_clusters_690) {
  points(combined_data[kmeans_result_selected_column_data_690$cluster == cluster_id, ]$XMax, 
         combined_data[kmeans_result_selected_column_data_690$cluster == cluster_id, ]$YMax, 
         col = "green", pch = 20, cex = .2)
  print(length(combined_data[kmeans_result_selected_column_data_690$cluster == cluster_id, ]$XMax))
  }

#570
plot(combined_data$XMax, 
     combined_data$YMax, 
     xlim = c(0, max(combined_data$XMax)), 
     ylim = c(0, max(combined_data$YMax)), 
     col = "lightgrey", 
     pch = 20, cex = .2, title(paste(titles[2])))
for (cluster_id in top_two_clusters_570) {
  points(combined_data[kmeans_result_selected_column_data_570$cluster == cluster_id, ]$XMax, 
         combined_data[kmeans_result_selected_column_data_570$cluster == cluster_id, ]$YMax, 
         col = "green", pch = 20, cex = .2)
  print(length(combined_data[kmeans_result_selected_column_data_570$cluster == cluster_id, ]$XMax))
  }

#520
plot(combined_data$XMax, 
     combined_data$YMax, 
     xlim = c(0, max(combined_data$XMax)), 
     ylim = c(0, max(combined_data$YMax)), 
     col = "lightgrey", 
     pch = 20, cex = .2, title(paste(titles[1])))
for (cluster_id in top_two_clusters_520) {
  points(combined_data[kmeans_result_selected_column_data_520$cluster == cluster_id, ]$XMax, 
         combined_data[kmeans_result_selected_column_data_520$cluster == cluster_id, ]$YMax, 
         col = "green", pch = 20, cex = .2)
  print(length(combined_data[kmeans_result_selected_column_data_520$cluster == cluster_id, ]$XMax))
}
