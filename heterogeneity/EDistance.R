### using E-distance for assessing heterogeneity for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20231220
### Last edited by Madeline E Melzer on 20240222

library(tidyverse)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(glue)
library(umap)
library(scperturbR)
library(vegan)
library(SingleCellExperiment)
library(ggpubr)
library(permute)
options(future.globals.maxSize = 4000 * 1024^2)

set.seed = 23

srcDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/"
eDistDataDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/data/"

datasets = c("FM01", 
             "FM02", 
             "FM03", 
             "FM04", 
             "FM05", 
             "FM06", 
             "FM08", 
             "Biorxiv",
             "non_cancer", 
             "ClonMapper", 
             "SPLINTR", 
             "LARRY",
             "cellTag",
             "watermelon",
             "TREX",
             "smartseq3_reads",
             "smartseq3_umis")

##########################################################################################################################
# calculating E-distance
##########################################################################################################################
# adapted from Goyal et al. 2023, E-dist.R

getEDistances <- function (seuratObject) {
  eDistances = data.frame(resolution = numeric(),
                          numClusters = integer(),
                          clusterA = character(),
                          clusterB = character(),
                          E_distance = numeric())
  
  eDistances_control = data.frame(resolution = numeric(),
                                  numClusters = integer(),
                                  clusterA = character(),
                                  clusterB = character(),
                                  E_distance = numeric())
  
  for (i in c(0.2, 0.4, 0.6, 0.8)) {
    resolution = i
    estats <- edist(seuratObject, groupby = paste0("clusters_res_", as.character(resolution))) #get EDistances
    estats_df <- as.data.frame(estats)
    
    # Set the row names as a column in estats_df
    estats_df$clusterA <- rownames(estats_df)
    
    # Pivot the data frame to a longer format
    estats_long <- estats_df %>%
      pivot_longer(-clusterA, names_to = "clusterB", values_to = "E_distance")
    
    # Add resolution and numClusters columns
    estats_long$resolution <- resolution
    estats_long$numClusters <- nrow(estats_df)
    
    # Append to eDistances
    eDistances <- bind_rows(eDistances, estats_long)
    
    ### Shuffling cluster numbers as a control
    
    cluster_data <- seuratObject[[paste0("clusters_res_", as.character(resolution))]]
    cluster_data_shuff = cluster_data
    cluster_data_shuff[[1]] <- sample(cluster_data[[1]])
    seuratObject@meta.data[, paste0("clusters_res_", as.character(resolution), "_shuffled")] <- cluster_data_shuff
    
    # Calculate E-distances for shuffled data
    estats_control <- edist(seuratObject, groupby = paste0("clusters_res_", as.character(resolution), "_shuffled"))
    estats_control_df <- as.data.frame(estats_control)
    estats_control_df$clusterA <- rownames(estats_control_df)
    estats_control_long <- estats_control_df %>%
      pivot_longer(-clusterA, names_to = "clusterB", values_to = "E_distance")
    estats_control_long$resolution <- resolution
    estats_control_long$numClusters <- nrow(estats_control_df)
    eDistances_control <- bind_rows(eDistances_control, estats_control_long)
    
  }
  print(eDistances)
  print(eDistances_control)
  return(list(eDistances = eDistances, eDistances_control = eDistances_control))
}


eDistResultsDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/results/"
eDistDataDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/data/"
singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/subsampledSingletObjects"


singletFiles <- list.files(singletObjectsDirectory, pattern = "\\.rds$", full.names = TRUE)

# Loop through each file and process
for(filepath in singletFiles) {
  # Extract dataset and sample from the file name
  parts <- strsplit(basename(filepath), "__")[[1]]
  dataset <- parts[1]
  sample <- gsub("_subsampledSinglets\\.rds$", "", parts[2])
  
  # Load the Seurat object
  seuratObject <- readRDS(filepath)
  
  # Run the eDistances function
  eDistances <- getEDistances(seuratObject)
  
  # Save the results
  write_csv(eDistances$eDistances, glue("{eDistDataDir}/{dataset}__{sample}__eDistance.csv"))
  write_csv(eDistances$eDistances_control, glue("{eDistDataDir}/{dataset}__{sample}__eDistance_control.csv"))
}

### for testing

seuratObject = readRDS("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/subsampledSingletObjects/Biorxiv__merged__subsampledSinglets.rds")
numNeighbors = 




#write_csv(eDistanceData, glue("{eDistResultsDir}/{dataset}__{sample}__eDistance.csv"))
#write_csv(eDistances_control, glue("{eDistResultsDir}/{dataset}__{sample}__eDistance_control.csv"))

##########################################################################################################################
# Aggregating all E-distance data
##########################################################################################################################


files <- list.files(eDistDataDir, pattern = "\\.csv$", full.names = TRUE)
processed_data <- tibble()

# Loop through each file
for (file in files) {
  parts <- strsplit(basename(file), "__")[[1]]
  dataset <- parts[1]
  sample <- gsub("(_eDistance.csv$|_eDistance_control.csv$)", "", parts[2])
  
  # Determine if the file is a control or sample based on the filename
  identity <- ifelse(grepl("control", file), "control", "sample")
  
  # Read the E-distance data
  eDistances_df <- read_csv(file) %>%
    mutate(clusterA = as.numeric(clusterA),
           clusterB = as.numeric(clusterB),
           identity = identity) %>%
    filter(!is.na(clusterA), !is.na(clusterB), clusterA != clusterB)
  
  # Add dataset and sample information
  eDistances_df <- mutate(eDistances_df, dataset = dataset, sample = sample)
  
  # Combine with the processed data
  processed_data <- bind_rows(processed_data, eDistances_df)
}

#write_csv(processed_data, glue("{eDistResultsDir}EDistances_allDataAndControls.csv"))

















avg_distances <- processed_data %>%
  group_by(dataset, numClusters) %>%
  #filter(numClusters == 4) %>%
  summarise(average_distance = mean(E_distance, na.rm = TRUE))

ggplot(avg_distances, aes(x = dataset, y = average_distance, fill = as.factor(numClusters))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Dataset", y = "Average Distance", fill = "Number of Clusters") +
  ggtitle("Average Distance vs Dataset for Each Number of Clusters")

plot <- ggplot(avg_distances, aes(x = dataset, y = average_distance, fill = as.factor(numClusters))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Dataset", y = "Average Distance", fill = "Number of Clusters") +
  ggtitle("Average Distance vs Dataset for Each Number of Clusters") +
  facet_wrap(~ numClusters, ncol = 4)  # Adjust ncol as needed to control the number of columns in the grid

# Print the plot
print(plot)










































##########################################################################################################################
# OLD
##########################################################################################################################
singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/"
eDistResultsDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/data"
eDistCombinedDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/results/"

singletCountsPath <- glue("{singletObjectsDirectory}singletCountBySample.csv")
singletCounts <- read_csv(singletCountsPath)

qualifiedPairs <- singletCounts %>% filter(numSinglets > 1000)

qualifiedPairs_splintr = qualifiedPairs %>% filter(dataset %in% c("TREX", "watermelon"))

for (i in 1:nrow(qualifiedPairs_splintr)) {
  set.seed(23)
  dataset <- qualifiedPairs_splintr$dataset[i]
  sample <- qualifiedPairs_splintr$sample[i]
  
  # Construct the path to the singlet RDS file
  singletFilePath <- glue("{singletObjectsDirectory}{dataset}/{sample}_singletsOnly.rds")
  
  # Check if the file exists before trying to read it
  if (file.exists(singletFilePath)) {
    object <- readRDS(singletFilePath)
    
    object <- getSubsampledSeurat(object, 1000)
    object <- getCellNeighborhood(object, 100)
    set.seed(45)
    object <- getCellNeighborhood(object, 100)
    eDistanceData <- getEDistances(object, 100)
    
    write_csv(eDistanceData, glue("{eDistResultsDir}/{dataset}__{sample}__eDistance.csv"))
  } else {
    message(glue("File {singletFilePath} does not exist."))
  }
}

# Assuming you have a list or vector of resolutions
resolutions <- c(0.2, 0.4, 0.6, 0.8)

# Initialize an empty dataframe to store the results
comparison_results <- tibble(dataset = character(), 
                             sample = character(), 
                             resolution = numeric(), 
                             cluster_pair_1 = character(), 
                             cluster_pair_2 = character(), 
                             distance_1 = numeric(), 
                             distance_2 = numeric())

comparison_results <- tibble()

for (i in 1:nrow(qualifiedPairs)) {
  dataset <- qualifiedPairs$dataset[i]
  sample <- qualifiedPairs$sample[i]
  
  eDistancesPath <- glue("{eDistResultsDir}/{dataset}__{sample}__eDistance.csv")
  
  # Check if the eDistances file exists
  if (file.exists(eDistancesPath)) {
    eDistances_df <- read_csv(eDistancesPath) %>%
      mutate(clusterA = as.numeric(clusterA),
             clusterB = as.numeric(clusterB)) %>%
      filter(!is.na(clusterA), !is.na(clusterB), clusterA != clusterB)
    
    # Loop through each resolution
    for (resolution in resolutions) {
      # Filter the eDistances dataframe for the current resolution
      eDistances_res_df <- filter(eDistances_df, resolution == resolution)
      
      # Randomly select two unique rows from the filtered dataframe for the current resolution
      if (nrow(eDistances_res_df) >= 3) {  # Ensure there are at least two rows to choose from
        selected_rows <- sample_n(eDistances_res_df, 3)
        
        # Add the dataset, sample, and resolution information to the selected rows
        selected_rows <- mutate(selected_rows, dataset = dataset, sample = sample)
        
        # Add the selected rows to the comparison_results dataframe
        comparison_results <- bind_rows(comparison_results, selected_rows)
      }
    }
  } else {
    message(glue("eDistances file {eDistancesPath} does not exist."))
  }
}

# View the final comparison results
print(comparison_results)
write_csv(comparison_results, glue("{eDistCombinedDir}threeSubsampledClustersPairwiseDistances.csv"))


avg_distances <- comparison_results %>%
  group_by(dataset, numClusters) %>%
  summarise(average_distance = mean(E_distance, na.rm = TRUE))



ggplot(avg_distances, aes(x = dataset, y = average_distance, fill = as.factor(numClusters))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  labs(x = "Dataset", y = "Average Distance", fill = "Number of Clusters") +
  ggtitle("Average Distance vs Dataset for Each Number of Clusters")


