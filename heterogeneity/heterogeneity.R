### assessing heterogeneity for Zhang Melzer et al 2024
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
#             MAIN             #
#                              #
################################
singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/"
#eDistResultsDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/data"
#eDistCombinedDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/results/"

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
  #object[["percent.mt"]] = PercentageFeatureSet(object, pattern = "^MT-")
  #object = subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
  
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
  pcaCoordinates = read_csv(glue("{singletObjectsDirectory}/pcaCoordinates/{dataset}__{sample}__pcaCoordinates.csv"))
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

#write_csv(shannonEntropy, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/shannonEntropy/shannonEntropy.csv")





################################
#                              #
#       Phenotypic Volume      #
#                              #
################################

# adapted from Fennell, Vassiliadis et al, 2022, scRNA_transcriptional_diversity_analysis.Rmd

singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/"
dataset = "FM01"
sample = "sample1"
object = readRDS(glue("{singletObjectsDirectory}{dataset}/{sample}_singletsOnly.rds"))

objsce <- as.SingleCellExperiment(object)
exp <- as.matrix(assays(objsce)$logcounts) #same as the data slot in the seurat object
meta <- colData(objsce)

testCovgene <- function(exp, meta) {
  #sub_meta=subset(meta, result.3 == cond )
  dat=exp[ , match(rownames(meta), colnames(exp))]
  dat=t(dat) #feature on col, cell on row
  dat=scale(dat) #standardize data scale by col/feature
  cols_with_na <- colSums(is.na(dat)) > 0
  dat_clean <- dat[, !cols_with_na]
  dat = dat_clean
  
  det=list()
  for (i in 1:50){ #permutation 
    set.seed(i)
    
    #cellsample=dat[1:nrow(dat), ]
    cellsample=sample(1:nrow(dat),1000, replace=F)
    
    gene_ranges <- apply(dat, 2, range)
    variable_genes <- apply(gene_ranges, 2, function(x) diff(x) != 0) # Identify genes with variation (where max != min)
    dat_variable <- dat[, variable_genes]
    
    genesample=sample(1:ncol(dat_variable), 1000, replace=F) #take long time to run so could sample genes, I also tried with all genes and it shows similar trend but bigger difference 
    
    cdat = cov(dat[cellsample,genesample]) #calculate gene-gene covariance matrix using sampled genes
    #cdat = cov(dat[cellsample,genesample]) #calculate gene-gene covariance matrix using all genes 
    
    svddat=svd(cdat) #singular value decomposition 
    
    eigen=svddat$d #get the eigen value 
    
    eigen=eigen[eigen!=0]
    
    logV=sum(log(eigen)) #log sum of non-zero eigen values
    #    logV=logV/length(genesample) #you could also norm by total number of genes 
    print(logV)
    det[[i]]=logV
  }
  
  return(det)
}

phenoVol = unlist(testCovgene(exp, meta))





















