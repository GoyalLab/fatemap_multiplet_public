### plotting clustering stability in doublet-containing vs not datasets for Zhang, Melzer et al 2024
### Created by Madeline E Melzer on 20240213
### Last edited by Madeline E Melzer on 20240311

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(ggpubr)
library(cowplot)
library(ggpattern)

set.seed(23)

barcodedDatasets = c("FM01", 
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
             "TREX")

barcodedRename = c("FM01" = "Goyal et al. 1", 
                   "FM02" = "Goyal et al. 2", 
                   "FM03" = "Goyal et al. 3", 
                   "FM04" = "Goyal et al. 4", 
                   "FM05" = "Goyal et al. 5", 
                   "FM06" = "Goyal et al. 6", 
                   "FM08" = "Goyal et al. 8", 
                   "non_cancer" = "Jiang et al.", 
                   #"ClonMapper", 
                   "Biorxiv" = "Jain et al.", 
                   #"SPLINTR", 
                   #"TREX",
                   #"LARRY",
                   "cellTag" = "CellTag-Multi",
                   "watermelon" = "Watermelon")

####################################################################################
# cluster 
####################################################################################

dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/clusteringStability/results"
plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/clusteringStability/"

files <- list.files(path = dataDirectory, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(files, read_csv)

# Combine all data frames into one
combined_data <- bind_rows(data_list)

# Find the number of clusters when there are no doublets for each dataset and resolution
baseline_clusters <- combined_data %>%
  filter(dbl_act == 0) %>%
  select(dataset, resolution, numClusters) %>%
  distinct()

# Merge this info back with the main data to compare with when there are doublets
comparison_data <- combined_data %>%
  left_join(baseline_clusters, by = c("dataset", "resolution"), suffix = c("", "_baseline"))

# Determine if the number of clusters is correct for the various doublet percentages
comparison_data <- comparison_data %>%
  mutate(correct_clusters = numClusters == numClusters_baseline)

#write_csv(comparison_data, paste0(dataDirectory, "clusteringStability_Louvain_allDatasets.csv"))
comparison_data = read.csv(paste0(dataDirectory, "clusteringStability_Louvain_allDatasets.csv"))

# Aggregate data for plotting
plot_data <- comparison_data %>%
  group_by(dataset, dbl_act, resolution, correct_clusters) %>%
  summarise()

plot_data$dataset <- factor(plot_data$dataset, levels = barcodedDatasets)

plot_data_filtered <- plot_data %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename)) %>%
  filter(resolution == 0.6)
  

# "heat map" of correct singlet identification
plot = ggplot(plot_data_filtered, aes(x = dataset, y = factor(dbl_act), fill = correct_clusters)) +
  geom_tile() +  # Use geom_tile to create the grid
  scale_fill_manual(values = c("TRUE" = "#7db0ea", "FALSE" = "#f6f2ee")) +  # Replace color1 and color2 with your color choices
  #facet_wrap(~ resolution, scales = "free_x", ncol = 1) +  # One column for each resolution
  labs(x = "", y = "Actual Doublet Rate (%)", fill = "Correct Number of Clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x labels for clarity
        strip.text.x = element_text(face = "bold"),
        legend.position = "bottom",
        axis.ticks = element_blank()) 
plot

ggsave(plot, file = paste0(plotDirectory, 'correctClustersHeatMap_0.6.svg'), width = 6, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'correctClustersHeatMap_0.6.png'), width = 6, height = 3)




plot_data_refiltered = plot_data_filtered %>% filter(correct_clusters == "TRUE")



########### plotting number of times clusters were correctly identifed vs doublet rate


correct_counts <- plot_data_filtered %>%
  filter(correct_clusters == TRUE) %>%  # Focus on correct identifications
  group_by(dbl_act) %>%
  summarise(correct_count = n())  # Count the number of correct identifications per dbl_act

# Bar graph
plot = ggplot(correct_counts, aes(x = dbl_act, y = correct_count)) +
  geom_bar(stat = "identity", fill = "#8cc53f") +
  #labs(x = "Doublet Rate", y = "Count of Correct Cluster Identifications", title = "Correct Cluster Identifications vs. Doublet Rate") +
  theme_classic()
plot

ggsave(plot, file = paste0(plotDirectory, 'barGraphOfNumCorrectClusterIdentifications_0.6.svg'), width = 2, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'barGraphOfNumCorrectClusterIdentifications_0.6.png'), width = 2, height = 4)



























########### plotting percent doublets per cluster vs dataset
comparison_data_subset = comparison_data %>% filter(dbl_act != 0, dbl_act == 10)

comparison_data_subset$resolution <- factor(comparison_data_subset$resolution, levels = c("0.2", "0.4", "0.6", "0.8"))
comparison_data_subset$dataset <- as.factor(comparison_data_subset$dataset)

# Create the boxplot
ggplot(comparison_data_subset, aes(x = dataset, y = pct_doublets, fill = resolution)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +  # Dodge positions the boxplots side-by-side within each dataset group
  scale_fill_brewer(palette = "Set1") +  # Use a predefined color palette for resolutions; adjust as needed
  labs(x = "Dataset", y = "Percentage of Doublets", fill = "Resolution") +  # Label axes and legend
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



#coloring the outliine by whether or not the true number of doublets is ther


ggplot(comparison_data_subset, aes(x = dataset, y = pct_doublets, fill = resolution)) +
  geom_boxplot(aes(color = correct_clusters), position = position_dodge(width = 0.75)) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_manual(values = c("TRUE" = "gold", "FALSE" = "black")) +  # 'NA' removes the outline for FALSE
  labs(x = "Dataset", y = "Percentage of Doublets", fill = "Resolution") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#### plotting dbl_act instead of for resolution

comparison_data = read.csv(paste0(dataDirectory, "clusteringStability_Louvain_allDatasets.csv"))

comparison_data_plot = comparison_data %>% filter(dbl_act != 0, resolution == 0.6)
comparison_data_plot$dbl_act <- as.numeric(as.character(comparison_data_plot$dbl_act))
comparison_data_plot$dbl_act <- factor(comparison_data_plot$dbl_act, levels = c(10, 20, 40))

# Create the boxplot
ggplot(comparison_data_plot, aes(x = dataset, y = pct_doublets, fill = dbl_act)) +
  geom_boxplot(aes(color = correct_clusters), position = 'dodge') +  # Use position_dodge to separate the boxplots for different dbl_act values
  scale_fill_brewer(palette = "Spectral") +  # Choose a color palette that suits your data; "Spectral" is just an example
  scale_color_manual(values = c("TRUE" = "#7db0ea", "FALSE" = "black")) +  # 'NA' removes the outline for FALSE
  labs(x = "Dataset", y = "Percentage of Doublets", fill = "Actual Doublet Rate (%)") +  # Label axes and legend
  theme_minimal() +  # Use a minimal theme for the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability





comparison_data$dataset <- as.factor(comparison_data$dataset)

# Filter the data
comparison_data_plot <- comparison_data %>% filter(dbl_act != 0, resolution == 0.6)

# Convert dbl_act to a numeric type first if it's not already
comparison_data_plot$dbl_act <- as.numeric(as.character(comparison_data_plot$dbl_act))

# Now convert dbl_act to a factor with the desired level order
comparison_data_plot$dbl_act <- factor(comparison_data_plot$dbl_act, levels = c(10, 20, 40))


comparison_data_plot$dataset <- factor(comparison_data_plot$dataset, levels = barcodedDatasets)
comparison_data_plot = comparison_data_plot %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename)) 

# Create the boxplot
plot = ggplot(comparison_data_plot, aes(x = dataset, y = pct_doublets, fill = dbl_act)) +
  geom_boxplot(position = position_dodge(width = 0.75), shape = 16) +
  scale_fill_brewer(palette = "Spectral") +
  scale_color_manual(values = c("TRUE" = "#7db0ea", "FALSE" = "black")) +
  labs(x = "Dataset", y = "Percentage of Doublets", fill = "Actual Doublet Rate (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
plot

ggsave(plot, file = paste0(plotDirectory, 'percentDoubletsPerCluster_res0.6.svg'), width = 7, height = 5)
ggsave(plot, file = paste0(plotDirectory, 'percentDoubletsPerCluster_res0.6.png'), width = 7, height = 5)



####### plotting UMAPS to figure out why this looks so odd

control = readRDS("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/functional_analysis_dataset/FM06/data/control.rds")
control = JoinLayers(control) # there are many samples which all need to be clustered together. 
### normalization
control = NormalizeData(control)
control = FindVariableFeatures(control)
control = ScaleData(control)

### dimension reduction
control <- RunPCA(object =control, verbose = FALSE)
control <- FindNeighbors(object=control, dims=1:50, reduction = "pca")
control <- FindClusters(object=control, resolution=0.6, algorithm = 1, verbose = FALSE)
control <- RunUMAP(object = control, reduction = "pca", dims = 1:50, verbose = FALSE)
DimPlot(object = control, reduction = "pca", group.by = "label")
DimPlot(object = control, reduction = "pca", group.by = "RNA_snn_res.0.6")



ten = readRDS("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/functional_analysis_dataset/Biorxiv/data/ten.rds")
ten = JoinLayers(ten) # there are many samples which all need to be clustered together. 
### normalization
ten = NormalizeData(ten)
ten = FindVariableFeatures(ten)
ten = ScaleData(ten)

### dimension reduction
ten <- RunPCA(object = ten, verbose = FALSE)
ten <- FindNeighbors(object = ten, dims = 1:50, reduction = "pca")
ten <- FindClusters(object = ten, resolution = 0.6, algorithm = 1, verbose = FALSE)
ten <- RunUMAP(object = ten, reduction = "pca", dims = 1:50, verbose = FALSE)
DimPlot(object = ten, reduction = "umap", group.by = "label")
DimPlot(object = ten, reduction = "umap", group.by = "RNA_snn_res.0.6")



twenty = readRDS("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/functional_analysis_dataset/FM06/data/twenty.rds")
twenty = JoinLayers(twenty) # there are many samples which all need to be clustered together. 
### normalization
twenty = NormalizeData(twenty)
twenty = FindVariableFeatures(twenty)
twenty = ScaleData(twenty)

### dimension reduction
twenty <- RunPCA(object = twenty, verbose = FALSE)
twenty <- FindNeighbors(object = twenty, dims = 1:50, reduction = "pca")
twenty <- FindClusters(object = twenty, resolution = 0.6, algorithm = 1, verbose = FALSE)
twenty <- RunUMAP(object = twenty, reduction = "pca", dims = 1:50, verbose = FALSE)
DimPlot(object = twenty, reduction = "pca", group.by = "label")
DimPlot(object = twenty, reduction = "pca", group.by = "RNA_snn_res.0.6")




forty = readRDS("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/functional_analysis_dataset/FM06/data/forty.rds")
forty = JoinLayers(forty) # there are many samples which all need to be clustered together. 
### normalization
forty = NormalizeData(forty)
forty = FindVariableFeatures(forty)
forty = ScaleData(forty)

### dimension reduction
forty <- RunPCA(object = forty, verbose = FALSE)
forty <- FindNeighbors(object = forty, dims = 1:50, reduction = "pca")
forty <- FindClusters(object = forty, resolution = 0.6, algorithm = 1, verbose = FALSE)
forty <- RunUMAP(object = forty, reduction = "pca", dims = 1:50, verbose = FALSE)
DimPlot(object = forty, reduction = "pca", group.by = "label")
DimPlot(object = forty, reduction = "pca", group.by = "RNA_snn_res.0.6")



