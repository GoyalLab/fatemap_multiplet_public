### determining doublet effects on clustering stability for Zhang Melzer et al. 2023
### Adapted from Xi and Li, 2021 Cell Systems
### Created by Madeline E Melzer on 20240103
### Last edited by Madeline E Melzer on 20240210

library(tidyverse)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(glue)
library(umap)
library(vegan)
library(SingleCellExperiment)
library(ggpubr)
options(future.globals.maxSize = 4000 * 1024^2)

set.seed(23)
###USED LOUVAIN ALGORITHM
####CALCULATE NUMBER OF CLUSTERS
#### CALCULATE SINGLET DISTRIBUTION PER CLUSTER IF THE REMAINING DROPLETS LED TO THE CORRECT NUMBER OF CELL CLUSTERS

getClusteredSeurat = function (seuratObject, numCells) {
  cluster_info <- data.frame(dataset = character(), sample = character(), dbl_pct = integer(), resolution = numeric(), num_clusters = integer())
  
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
  analysis.seurat.subset <- FindNeighbors(object=analysis.seurat.subset, dims=1:50, reduction = "pca")
  
  # save PCA coordinates
  pcaCoordinates = analysis.seurat.subset@reductions[['pca']]@cell.embeddings
  cells_PCA = rownames(pcaCoordinates) #CellIds with Sample number as prefix
  pcaCoordinates = as_tibble(pcaCoordinates) #UMAP coordinates in tibble
  pcaCoordinates = pcaCoordinates %>% mutate(cellID = cells_PCA)
  write_csv(pcaCoordinates, file= glue("{clusteringStabilityDirectory}pcaCoordinates/{dataset}__{condition}__pcaCoordinates__{dbl_act}.csv"))
  
  # cluster at different resolutions
  for (i in c(0.2, 0.4, 0.6, 0.8)) {
    resolution = i
    analysis.seurat.subset <- FindClusters(object=analysis.seurat.subset, resolution=resolution, algorithm = 1, verbose = FALSE) #clusters based on Louvain algorithm (algorithm = 1)
    cluster_column_name <- glue("SCT_snn_res.{resolution}")
    num_clusters <- length(unique(analysis.seurat.subset@meta.data[[cluster_column_name]]))
    analysis.seurat.subset <- RunUMAP(object = analysis.seurat.subset, reduction = "pca", dims = 1:50, verbose = FALSE)
    cluster_info <- rbind(cluster_info, data.frame(dataset = dataset, sample = sample, dbl_pct = 0.0, resolution = resolution, numClusters = num_clusters))
  }
  #write_csv(cluster_info, glue("{resultsDirectory}{dataset}__{sample}__cellClusters.csv"))
  analysis.seurat.subset
}


##########################################################################################################################
# proportion of singlets in correctly identified clusters (cluster quality) 
##########################################################################################################################
functionalObjectsDirectory = "/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset/"
clusteredObjectsDirectory = "/projects/p31666/melzer/ZhangMelzerEtAl/clusteringStability/clusteredObjects"
resultsDirectory = "/projects/p31666/melzer/ZhangMelzerEtAl/clusteringStability/results/"


#functionalObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/functional_analysis_dataset/"
#clusteredObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/clusteringStability/clusteredObjects"
#resultsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/clusteringStability/results/"

set.seed(23)
datasets <- list.dirs(functionalObjectsDirectory, full.names = FALSE, recursive = FALSE)
datasets  = c("FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "LARRY", "non_cancer", "smartseq3", "smartseq3_reads", "smartseq3_umis", "SPLINTR", "TREX", "watermelon")

for (dataset in datasets) {
  results <- data.frame(dataset = character(), dbl_act = numeric(), resolution = numeric(), numClusters = integer(), cluster = integer(), num_singlets = numeric(), num_doublets = numeric(), pct_doublets = numeric(), stringsAsFactors = FALSE)
  conditions <- c("control", "forty", "twenty", "ten")
  dbl_act_values <- c(0, 40, 20, 10)
  condition_path <- file.path(functionalObjectsDirectory, dataset, "data")
  
  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    dbl_act <- dbl_act_values[i]
    
    file_path <- file.path(condition_path, paste0(condition, ".rds"))
    
    if (file.exists(file_path)) {
      object <- readRDS(file_path)
      object <- getClusteredSeurat(object, 1000)  # Ensure this function clusters and adds cluster info to meta.data
      
      # Loop through resolutions
      for (resolution in c(0.2, 0.4, 0.6, 0.8)) {
        cluster_column_name <- paste0("SCT_snn_res.", resolution)
        num_clusters <- length(unique(object@meta.data[[cluster_column_name]]))
        saveRDS(object, glue("{clusteredObjectsDirectory}/{dataset}__{condition}__clustered.rds"))
        # Calculate percentage of doublets for each cluster
        for (cluster in unique(object@meta.data[[cluster_column_name]])) {
          cluster_cells <- object@meta.data[object@meta.data[[cluster_column_name]] == cluster, ]
          num_doublets <- sum(cluster_cells$label == "doublet", na.rm = TRUE)
          num_singlets <- sum(cluster_cells$label == "singlet", na.rm = TRUE)
          pct_doublets <- mean(cluster_cells$label == "doublet", na.rm = TRUE) * 100  # Assuming 'label' column exists
          
          # Append the information to the results dataframe
          results <- rbind(results, data.frame(dataset = dataset, 
                                               dbl_act = dbl_act, 
                                               resolution = resolution, 
                                               numClusters = num_clusters, 
                                               cluster = cluster, 
                                               num_doublets = num_doublets, 
                                               num_singlets = num_singlets, 
                                               pct_doublets = pct_doublets))
        }
      }
      
    } else {
      message(glue("File {file_path} does not exist."))
    }
  }
  write.csv(results, file = glue("{resultsDirectory}{dataset}__clusteringStability.csv"), row.names = FALSE)
}

# Write results to CSV
write.csv(results, file = glue("{resultsDirectory}allDatasetsClusteringStability.csv"), row.names = FALSE)



##########################################################################################################################
# making average plot
##########################################################################################################################



