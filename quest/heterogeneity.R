### assessing heterogeneity for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20231220
### Last edited by Madeline E Melzer on 20240210

library(tidyverse)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(glue)
library(umap)
#library(scperturbR)
library(vegan)
library(SingleCellExperiment)
library(ggpubr)
options(future.globals.maxSize = 4000 * 1024^2)

set.seed = 23

srcDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/"

barcoded = c("FM01", 
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

################################
#                              #
#  Making Singlet-Only Objects #
#                              #
################################

### function for extracting singlets
extractSinglets = function(dataset, sample) {
  rnaDirectory = glue("{srcDirectory}/{dataset}/10X/{sample}/")
  singletsListDirectory = glue("{srcDirectory}/{dataset}/fatemapID/")
  singlets_file = glue("{singletsListDirectory}{dataset}_singlets_all.txt")
  
  expression_matrix <- Read10X(data.dir = rnaDirectory)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  
  # Loop through the numbers -1 to -6 to check for each suffix
  for (i in 1:6) {
    suffix <- paste0("-", i)  # Create the suffix string ("-1", "-2", ..., "-6")
    
    # Check if the current suffix is part of the 10X barcodes
    if (grepl(suffix, seurat_object |> colnames() |> head(1))) {
      # Read singlets and append the current suffix
      singlets <- read.delim(singlets_file, header = FALSE, sep = "\n")[, 1] |>
        as.vector() |>
        paste0(suffix)
      break  # Exit the loop once a matching suffix is found
    }
  }
  
  # If no matching suffix was found, remove potential suffixes "-1" to "-6"
  if (is.null(singlets)) {
    singlets <- read.delim(singlets_file, header = FALSE, sep = "\n")[, 1] |>
      as.vector()
    for (i in 1:6) {
      singlets <- singlets |> str_remove(paste0("-", i))
    }
  }

  #subset Seurat object so it only contains singlets
  singletsSeurat = subset(seurat_object, cells = singlets)
  
  # Count the number of singlets present in the Seurat object
  numSinglets <- length(x = colnames(singletsSeurat))
  
  singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/"
  # Create dataset directory if it does not exist yet
  datasetSingletDirectory <- glue("{singletObjectsDirectory}{dataset}/")
  if (!dir.exists(datasetSingletDirectory)) {
    dir.create(datasetSingletDirectory, recursive = TRUE)
  }
  
  singletObjectPath<- glue("{datasetSingletDirectory}{sample}_singletsOnly.rds")
  # Save the RDS file if it does not exist yet
  if (!file.exists(singletObjectPath)) {
    saveRDS(singletsSeurat, file = singletObjectPath)
  }
  
  return(list(singletsSeurat = singletsSeurat, numSinglets = numSinglets))
}


### making all the singlet RDS objects and counting singlets per sample
singletCounts = data.frame(dataset = character(), sample = character(), numSinglets = integer(), stringsAsFactors = FALSE)

datasetDirs = c("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/Biorxiv/") #can be list of all dires, watermelon added later

# Loop through each dataset directory
for (datasetDir in datasetDirs) {
  datasetName = basename(datasetDir)
  rnaDir = file.path(datasetDir, "10X")
  
  if (dir.exists(rnaDir)) {
    sampleDirs = list.dirs(rnaDir, recursive = FALSE, full.names = TRUE)
    for (sampleDir in sampleDirs) {
      sampleName = basename(sampleDir)
      
      # Extract singlets for current dataset and sample
      singletCount = extractSinglets(datasetName, sampleName)
      
      # Append the dataset, sample, and number of singlets to the results data frame
      singletCounts = rbind(singletCounts, data.frame(dataset = datasetName, sample = sampleName, numSinglets = singletCount$numSinglets))
      
      print(glue("Dataset: {datasetName}, Sample: {sampleName}, Number of Singlets: {singletCount$numSinglets}"))
    }
  }
}

#write.csv(singletCounts, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/singletCountBySample.csv", row.names= FALSE)


### merge Jain et al. ("Biorxiv") DMSO A and B to make it have over 1000 cells (only 800 cells in smallest sample)

sample1_dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/Biorxiv/10X/1_DMSO_A"
sample2_dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/Biorxiv/10X/2_DMSO_B"

dataset = "Biorxiv"
sample = "1_DMSO_A"
#setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("//Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/fatemapID/")
#sample1 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.1)
sample1 = Read10X(data_dir)
sample1 = CreateSeuratObject(sample1)
sample1@meta.data$sample <- "1_DMSO_A"

sample = "2_DMSO_B"
#setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/fatemapID/")
#sample2 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.1)
sample2 = Read10X(data_dir)
sample2 = CreateSeuratObject(sample2)
sample2@meta.data$sample <- "2_DMSO_B"

#saveRDS(sample1, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/sample1/sample1.rds")
#saveRDS(sample2, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/sample2/sample2.rds")

sample1 = PercentageFeatureSet(sample1, pattern = "^MT-", col.name = "percent.mt")
sample2 = PercentageFeatureSet(sample2, pattern = "^MT-", col.name = "percent.mt")

s1s2_scTransform <- merge(sample1, y = sample2, add.cell.ids = c("S1", "S2"), project = "S1S2")

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "sample")
s1s2_scTransform.list <- s1s2_scTransform.list[c("1_DMSO_A", "2_DMSO_B")]

for (i in 1:length(s1s2_scTransform.list)) {
  s1s2_scTransform.list[[i]] <- SCTransform(s1s2_scTransform.list[[i]])
}

s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 7000)

s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, 
                                            verbose = FALSE)

s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1s2_scTransform.features, verbose = FALSE)
s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", verbose = FALSE)

s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

# saving
setwd("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/Biorxiv/10X/DMSO_AB/")
matrix <- s1s2_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1s2_scTransform.integrated@meta.data$label
samples = s1s2_scTransform.integrated@meta.data$sample

Matrix::writeMM(matrix, file = "matrix.mtx")
#system("gzip matrix.mtx")
write.table(features, file = gzfile("features.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile("barcodes.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)








################################
#                              #
#             MAIN             #
#                              #
################################
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




###

processed_data <- tibble()

# Loop through each qualified dataset-sample pair
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
    
    # Add dataset and sample information
    eDistances_df <- mutate(eDistances_df, dataset = dataset, sample = sample)
    
    # Combine with the processed data
    processed_data <- bind_rows(processed_data, eDistances_df)
  } else {
    message(glue("eDistances file {eDistancesPath} does not exist."))
  }
}

write_csv(processed_data, glue("{eDistCombinedDir}allClusterPairwiseDistances.csv"))


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


################################
#                              #
#       Clustering             #
#                              #
################################

singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/"
dataset = "FM01"
sample = "sample1"
object = readRDS(glue("{singletObjectsDirectory}{dataset}/{sample}_singletsOnly.rds"))

### filter, subsample, normalize, dimension reduction, and multi-res clustering of Seurat object 
## (clusters_res_{resolution} in meta.data)
## returns labeled seurat object and pcaCoordinates
getSubsampledSeurat <- function (seuratObject, numCells) {
  ### filter
  object[["percent.mt"]] = PercentageFeatureSet(object, pattern = "^MT-")
  object = subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
  
  ### subsample numCells cells from the sample
  analysis.seurat = object
  cells <- sample(colnames(analysis.seurat), size =numCells, replace=F)
  analysis.seurat.subset <- analysis.seurat[, cells]
  
  ### normalization
  analysis.seurat.subset <- SCTransform(analysis.seurat.subset, verbose = FALSE) # normalizing, scaling, finding variable features
  
  ### dimension reduction
  analysis.seurat.subset <- RunPCA(object = analysis.seurat.subset, assay = "SCT", verbose = FALSE)
  analysis.seurat.subset <- FindNeighbors(object=analysis.seurat.subset, dims=1:50, reduction = "pca", force.recalc = TRUE) #neighborhood overlap, jaccard index
  
  # save PCA coordinates
  pcaCoordinates = analysis.seurat.subset@reductions[['pca']]@cell.embeddings
  cells_PCA = rownames(pcaCoordinates) #CellIds with Sample number as prefix
  pcaCoordinates = as_tibble(pcaCoordinates) #UMAP coordinates in tibble
  pcaCoordinates = pcaCoordinates %>% mutate(cellID = cells_PCA)
  write_csv(pcaCoordinates, file= glue("{singletObjectsDirectory}/pcaCoordinates/{dataset}__{sample}__pcaCoordinates.csv"))
  
  # cluster at different resolutions
  for (i in c(0.2, 0.4, 0.6, 0.8)) {
    resolution = i
    print(as.character(resolution))
    cluster_col_name <- paste0("clusters_res_", as.character(resolution))
    analysis.seurat.subset <- FindClusters(object=analysis.seurat.subset, resolution=resolution, verbose = FALSE) #clusters based on their SNN
    colnames(analysis.seurat.subset@meta.data)[which(colnames(analysis.seurat.subset@meta.data) == "seurat_clusters")] <- cluster_col_name
    analysis.seurat.subset <- RunUMAP(object = analysis.seurat.subset, reduction = "pca", dims = 1:50, verbose = FALSE)
    #print(DimPlot(analysis.seurat.subset, reduction = "umap", group.by = cluster_col_name) + ggtitle(paste("UMAP at Resolution", resolution)))
  }
  
  analysis.seurat.subset
}

### subsample cells randomly (across all UMAP space) 
## (subsampled in meta.data)
## returns labeled Seurat object
getCellDiaspora <- function (seuratObject, numCells) {
  cells <- sample(colnames(seuratObject), size = numCells, replace=F)
  seuratObject$subsampled <- ifelse(colnames(seuratObject) %in% cells, "subsampled", "notSubsampled")
  
  print(DimPlot(seuratObject, reduction = "umap", group.by = "subsampled") +
    ggtitle("UMAP with Subsampled Cells Highlighted"))
  
  seuratObject
}

### subsample cells that neighbor a randomly chosen cell 
## (nn_label in meta.data)
## returns labeled Seurat object

getCellNeighborhood <- function (seuratObject, numNeighbors) {
  chosenCell <- sample(colnames(seuratObject), size = 1)
  pcaCoordinates = read_csv(glue("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/pcaCoordinates/{dataset}__{sample}__pcaCoordinates.csv"))
  chosenCell_coords <- pcaCoordinates %>%
    dplyr::filter(cellID == chosenCell) %>%
    dplyr::select(-cellID)
  
  # calculate distances from the chosen cell to all others
  distance_matrix <- as.matrix(dist(rbind(chosenCell_coords, pcaCoordinates[, -ncol(pcaCoordinates)])))
  distances <- distance_matrix[1, -1] # The first row contains distances from the chosen cell to all others
  
  # select top 200 nearest neighbors based on the smallest distances
  nearest_neighbors <- sort(distances, index.return = TRUE)$ix[1:numNeighbors]
  nn_cells <- pcaCoordinates$cellID[nearest_neighbors]
  
  # label these nearest neighbors in your Seurat object
  nn_label_name <- if ("nn_label" %in% names(seuratObject@meta.data)) {
    "nn_label2"
  } else {
    "nn_label"
  }
  seuratObject@meta.data[[nn_label_name]] <- ifelse(colnames(seuratObject) %in% nn_cells, "neighbor", "other")
  
  # visualize the results using UMAP and highlight the nearest neighbors
 print(DimPlot(seuratObject, reduction = "umap", group.by = "nn_label") +
    ggtitle("UMAP with Nearest Neighbors Highlighted") )
 print(DimPlot(seuratObject, reduction = "pca", group.by = "nn_label") +
         ggtitle("PCA with Nearest Neighbors Highlighted") )
  
  seuratObject
}

################################
#                              #
#       Shannon Entropy        #
#                              #
################################

# adapted from Goyal et al. 2023, FM01_ShannonEntropy_barcodes.R  and Fennell, Vassiliadis et al., 2022, scRNA_transcriptional_diversity_analysis.Rmd

###Functionalize Count normalized matrix
object = read_rds("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/FM01/sample1_singletsOnly.rds")
datset = "FM01"
sample = "sample1"
object = getSubsampledSeurat(object, 1000)

singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/"
resolutions <- c(0.2, 0.4, 0.6, 0.8)  # Update this list as needed
shannonEntropy <- data.frame()

qualifiedPairs_subset = qualifiedPairs %>% filter(dataset %in% c("SPLINTR", "TREX", "watermelon"))

for (i in 1:nrow(qualifiedPairs_subset)) {
  dataset <- qualifiedPairs_subset$dataset[i]
  sample <- qualifiedPairs_subset$sample[i]
  object = readRDS(glue("{singletObjectsDirectory}{dataset}/{sample}_singletsOnly.rds"))
  object = getSubsampledSeurat(object, 1000)
  object = getCellDiaspora(object, 200)
  saveRDS(object, glue("{singletObjectsDirectory}{dataset}/{sample}_singletsOnly_normSub.rds"))
  
  for (resolution in resolutions) {
    seurat_clusters <- paste0("clusters_res_", resolution)
    cluster_counts <- table(object@meta.data[[seurat_clusters]])
    cluster_matrix <- matrix(cluster_counts, ncol = 1)
    shannon_diversity_cluster <- diversity(cluster_matrix, index = "shannon")
    
    spread_counts = table(object@meta.data$subsampled)
    spread_matrix <- matrix(spread_counts, ncol = 1)
    shannon_diversity_spread <- diversity(spread_matrix, index = "shannon")
    
    temp_df <- data.frame(
      dataset = dataset,
      sample = sample,
      resolution = resolution,
      clusterDiversity = shannon_diversity_cluster,
      spreadDiversity = shannon_diversity_spread
    )
    
    shannonEntropy <- rbind(shannonEntropy, temp_df)
    
  }
}

write_csv(shannonEntropy, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/shannonEntropy/shannonEntropy.csv")

long_data <- pivot_longer(shannonEntropy, 
                          cols = c(clusterDiversity, spreadDiversity),
                          names_to = "DiversityType", 
                          values_to = "DiversityValue")

long_data_filtered <- long_data %>% filter(DiversityValue > 0.51)


# Plotting Shannon Diversity Index
plot <- ggplot(long_data_filtered, aes(y = DiversityValue, x = dataset, col = dataset)) + 
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, shape = 16) +
  stat_compare_means(size = 4, label.x.npc = "left", paired = TRUE) +
  facet_wrap(~resolution) +
  labs(x = "", y = "Shannon Diversity Index") +
  theme_classic() +
  theme(legend.position = "none")
plot

#Dataframe
#confirm order of barcodes is same/paired:
#all(entropy_random$BC50StarcodeD8==entropy_clones$BC50StarcodeD8)
plot <- ggplot(long_data_filtered, aes(y = DiversityValue, x = dataset, col = dataset)) + 
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, shape = 16) +
  stat_compare_means( size =4, label.x.npc = "left", paired = T)+
  facet_wrap(.~resolution)+
  labs(x="", y = "Shannon Diversity Index")+
  theme_classic() + 
  theme(legend.position="none") 
plot
ggsave(plot = plot,filename=paste0(plotDirectory, "ShannonDiversityIndex_resolutions.svg"), width= 18, height =12)

##############
#For normalized entropy - Shannon's Equitability Index
### I dont think I can do this because of the fact i do not have numClusters!!!!!!!! BAH
plot <-ggplot(long_data_filtered, aes(y = DiversityValue/log(nrow(DiversityValue)-1), x = dataset, col = dataset)) + 
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, shape = 16) +
  stat_compare_means( size =4, label.x.npc = "left", paired = T)+
  facet_wrap(.~resolution)+
  labs(x="", y = "Shannon Diversity Index")+
  theme_classic() + 
  theme(legend.position="none") 
plot

ggsave(plot = plot,filename=paste0(plotDirectory, "ShannonEquitabilityIndex.svg"), width= 12, height = 12)








################################
#                              #
#       Phenotypic Volume      #
#                              #
################################

# adapted from Fennell, Vassiliadis et al, 2022, scRNA_transcriptional_diversity_analysis.Rmd


testCovgene <- function(dat_variable, numCells) {
  #sub_meta=subset(meta, result.3 == cond )
  
  det=list()
  for (i in 1:10){ #permutation 
    set.seed(i)
    
    cellsample=sample(1:nrow(dat_variable), numCells, replace=F)
    
    cdat = cov(dat_variable[cellsample, ]) #calculate gene-gene covariance matrix using sampled genes
    
    svddat=svd(cdat) #singular value decomposition 
    
    eigen=svddat$d #get the eigen value 
    
    eigen=eigen[eigen!=0]
    
    logV=sum(log(eigen)) #log sum of non-zero eigen values
    print(logV)
    det[[i]]=logV
  }
  
  return(det)
}






pvResultsDir = "~/ZhangMelzerEtAl/data/phenotypicVolume/data_variableGenes"
pvCombinedDir = "~/ZhangMelzerEtAl/data/phenotypicVolume/results/"
singletObjectsDirectory = "/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/subsampledSingletObjects"

singletFiles <- list.files(singletObjectsDirectory, pattern = "\\.rds$", full.names = TRUE)

#biorxiv = readRDS("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/subsampledSingletObjects/Biorxiv__merged__subsampledSinglets.rds")
#DefaultAssay(biorxiv) <- 'RNA'

#FM01 = readRDS("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/subsampledSingletObjects/FM01__sample1__subsampledSinglets.rds")
#object = FM01
#numCells = 1000

singletFiles <- singletFiles[-1]  # Remove the first file from the list (this is Biorxiv)

#singletFiles = c("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/subsampledSingletObjects/Biorxiv__merged__subsampledSinglets.rds")

for(filepath in singletFiles) {
  # Extract dataset and sample from the file name
  print(filepath)
  parts <- strsplit(basename(filepath), "__")[[1]]
  dataset <- parts[1]
  sample <- gsub("_subsampledSinglets\\.rds$", "", parts[2])
  
  seuratObject <- readRDS(filepath)
  #DefaultAssay(seuratObject) <- 'RNA' # added for biorxiv merged dataset ONLY because is a result of scTransform of 2 samples
  #seuratObject = JoinLayers(seuratObject) # added for biorxiv merged dataset ONLY because is a result of scTransform of 2 samples
  #seuratObject = NormalizeData(seuratObject) # added for biorxiv merged dataset ONLY because is a result of scTransform of 2 samples
  #seuratObject = FindVariableFeatures(seuratObject) # added for biorxiv merged dataset ONLY because is a result of scTransform of 2 samples
  #seuratObject = ScaleData(seuratObject) # added for biorxiv merged dataset ONLY
  
  objsce <- as.SingleCellExperiment(seuratObject)
  exp <- as.matrix(assays(objsce)$logcounts) #same as the data slot in the seurat object
  meta <- colData(objsce)
  
  dat=exp[ , match(rownames(meta), colnames(exp))]
  dat=t(dat) #feature on col, cell on row
  #print(dat)
  dat = scale(dat)
  #standardize data scale by col/feature
  #cols_with_na <- colSums(is.na(dat)) > 0 # dont need to make data "clean" if we are using the 1000 most variable features
  #dat_clean <- dat[, !cols_with_na]
  #dat = dat_clean
  
  #selecting the 1000 most variable genes from the object to use for the following analysis
  variable_genes <- head(VariableFeatures(seuratObject), 1000)
  dat_variable <- dat[, variable_genes, drop = FALSE]
  
  allPV_values = unlist(testCovgene(seuratObject, 1000))
  
  seuratObject <- getCellNeighborhood(seuratObject, 100) #this is for the control
  object_subset <- subset(seuratObject, subset = nn_label == "neighbor")
  neighborhoodPV_values = unlist(testCovgene(object_subset, 100))
  
  pv_df <- data.frame(allPV = allPV_values, neighborhoodPV = neighborhoodPV_values)
  write_csv(pv_df, glue("{pvResultsDir}/{dataset}__{sample}__phenotypicVolume.csv"))
}














seuratObject <- RunPCA(object = seuratObject, verbose = FALSE)
seuratObject <- FindNeighbors(object=seuratObject, dims=1:50, reduction = "pca", force.recalc = TRUE)
seuratObject <- RunUMAP(object = seuratObject, reduction = "pca", dims = 1:50, verbose = FALSE)
print(DimPlot(seuratObject, reduction = "umap", group.by = "sample") + ggtitle(paste("UMAP at Resolution", resolution)))






singletCountsPath <- glue("{singletObjectsDirectory}singletCountBySample.csv")
singletCounts <- read_csv(singletCountsPath)
qualifiedPairs <- singletCounts %>% filter(numSinglets > 1000)

qualifiedPairs_filtered = qualifiedPairs %>% filter(dataset %in% c("SPLINTR", "TREX", "watermelon"))

for (i in 1:nrow(qualifiedPairs_filtered)) {
  set.seed(23)
  dataset <- qualifiedPairs_filtered$dataset[i]
  sample <- qualifiedPairs_filtered$sample[i]
  
  # Construct the path to the singlet RDS file
  singletFilePath <- glue("{singletObjectsDirectory}{dataset}/{sample}_singletsOnly.rds")
  
  # Check if the file exists before trying to read it
  if (file.exists(singletFilePath)) {
    object <- readRDS(singletFilePath)
    object <- getSubsampledSeurat(object, 1000)
    allPV_values = unlist(testCovgene(object, 1000))
    object <- getCellNeighborhood(object, 100)
    object_subset <- subset(object, subset = nn_label == "neighbor")
    neighborhoodPV_values = unlist(testCovgene(object_subset, 100))
    
    pv_df <- data.frame(allPV = allPV_values, neighborhoodPV = neighborhoodPV_values)
    write_csv(pv_df, glue("{pvResultsDir}/{dataset}__{sample}__phenotypicVolume.csv"))
    
  } else {
    message(glue("File {singletFilePath} does not exist."))
  }
}







#plot
library(ggplot2)
library(ggpubr)
library(cowplot)

plotdata<-data.frame(volume=phenoVol,
                     cell=c(rep(c("phenoVol"), times=c(length(phenoVol)))))
head(plotdata)
dim(plotdata)

ggplot(plotdata, aes(x=cell, y=volume, color=cell))+geom_violin()+geom_boxplot(width=0.1)+
  #        stat_compare_means()+ #this would give mann whitney test result
  theme_cowplot(12)+ylab("log of phenotypic volume(gene-gene)")+
  theme(legend.position="none") + 
  ggpubr::stat_compare_means(ref.group = "comb", method = "t.test", method.args = list(alternative = "less"))
ggsave("kras_v_comb_PV_genev2.pdf")





















################################
#                              #
#       E-distance             #
#                              #
################################

# adapted from Goyal et al. 2023, E-dist.R
eDistResultsDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/data/"
singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/"
dataset = "SPLINTR"
sample = "chemoDay2_1"
splintr = readRDS(glue("{singletObjectsDirectory}{dataset}/{sample}_singletsOnly.rds"))

splintr = getSubsampledSeurat(splintr, 1000)
splintr = getCellNeighborhood(splintr, 100)
splintr = getCellNeighborhood(splintr, 100)
seuratObject = splintr
eDistanceData = getEDistances(splintr, 100)

write_csv(eDistanceData, glue("{eDistResultsDir}/{dataset}__{sample}__eDistance.csv"))

getEDistances <- function (seuratObject, numNeighbors) {
  eDistances = data.frame(resolution = numeric(),
                          numClusters = integer(),
                          clusterA = character(),
                          clusterB = character(),
                          E_distance = numeric())

  for (i in c(0.2, 0.4, 0.6, 0.8)) {
    resolution = i
    estats <- edist(seuratObject, groupby = paste0("clusters_res_", as.character(resolution)))
    print(estats)
    
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
  }
  
  seuratObject$nn_label_cross <- apply(seuratObject@meta.data[, c("nn_label", "nn_label2")], 1, function(x) {
    if (x[1] == "neighbor" & x[2] == "neighbor") {
      return("Neither")
    } else if (x[1] == "neighbor") {
      return("Label1_Only")
    } else if (x[2] == "neighbor") {
      return("Label2_Only")
    } else {
      return("Neither")
    }
  })
  
  edist_nn_label_cross <- edist(seuratObject, groupby = "nn_label_cross", reduction = "pca")
  estats_df <- as.data.frame(edist_nn_label_cross)
  
  # Set the row names as a column in estats_df
  estats_df$clusterA <- rownames(estats_df)
  
  # Pivot the data frame to a longer format
  estats_long <- estats_df %>%
    pivot_longer(-clusterA, names_to = "clusterB", values_to = "E_distance")
  
  # Add resolution and numClusters columns
  estats_long$resolution <- numNeighbors
  estats_long$numClusters <- nrow(estats_df)
  
  # Append to eDistances
  eDistances <- bind_rows(eDistances, estats_long)
  
  eDistances
}


