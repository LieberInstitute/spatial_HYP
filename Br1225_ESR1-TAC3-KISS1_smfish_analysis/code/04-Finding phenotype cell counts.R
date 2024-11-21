# Load necessary library
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
# Set working directory to where the last files were generated from
setwd("/Users/atharv.chandra/Desktop/HYP_1225/Two Cluster Output")
pre_final_analysis_path <- "/Users/atharv.chandra/Desktop/HYP_1225/Analysis"
copy_count_name <- "HYP1225_selected_Three_Clusters.xlsx"

# List all .xlsx files in the folder
file_list <- list.files(pattern = "\\.xlsx$", full.names = TRUE)

# Santity check
file_list 

# Function to read the "phenotype" column from an Excel file and count occurrences
read_and_count_phenotypes <- function(file_path) {
  # Read the Excel file
  data <- read_excel(file_path)
  
  # Check if "phenotype" column exists
  if ("phenotype" %in% colnames(data)) {
    # Count occurrences of each phenotype
    count_data <- data %>%
      count(phenotype, name = "count") %>%
      mutate(file_name = basename(file_path))  # Add file name as a column
    return(count_data)
  } else {
    warning(paste("Column 'phenotype' not found in file:", file_path))
    return(NULL)
  }
}

# Read and count phenotypes from each file
phenotype_counts_list <- lapply(file_list, read_and_count_phenotypes)

# Combine all counts into a single data frame
combined_counts <- bind_rows(phenotype_counts_list)

# Reshape data to have phenotypes as rows and file names as columns
phenotype_summary <- combined_counts %>%
  pivot_wider(names_from = file_name, values_from = count, values_fill = list(count = 0)) %>%
  arrange(phenotype)

# Display the resulting data frame
print(phenotype_summary)

# Save the count data in a new excel file
write_xlsx(phenotype_summary, path = file.path(pre_final_analysis_path, paste("04-",copy_count_name)))

