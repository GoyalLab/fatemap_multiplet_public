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
library(dbscan)
options(future.globals.maxSize = 4000 * 1024^2)

set.seed(23)

##########################################################################################################################
# proportion of singlets in correctly identified clusters (cluster quality) 
##########################################################################################################################
functionalObjectsDirectory = "/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset"
#functionalObjectsDirectory = "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow" #OVERFLOW OBJECTS
clusteredObjectsDirectory = "/projects/p31666/melzer/ZhangMelzerEtAl/clusteringStability/objects"
pcaCoordinatesDirectory = "/projects/p31666/melzer/ZhangMelzerEtAl/clusteringStability/pcaCoordinates"
resultsDirectory = "/projects/p31666/melzer/ZhangMelzerEtAl/clusteringStability/results"

#functionalObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/functional_analysis_dataset/"
#clusteredObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/clusteringStability/clusteredObjects"
#resultsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/clusteringStability/results/"

set.seed(23)
#datasets <- list.dirs(functionalObjectsDirectory, full.names = FALSE, recursive = FALSE)

datasets_initial = c("Biorxiv", "FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "non_cancer","smartseq3_umis", "watermelon")
datasets_overflow  = c("cellTag", "ClonMapper","LARRY", "smartseq3_reads", "SPLINTR", "TREX")

datasets_initial = c("Biorxiv")
dbl_act = 40
condition = "forty"

for (dataset in datasets_initial) {
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
      object = JoinLayers(object) # there are many samples which all need to be clustered together. 
      
      ### normalization
      object = NormalizeData(object)
      object = FindVariableFeatures(object)
      object = ScaleData(object)
      
      ### dimension reduction
      object <- RunPCA(object = object, verbose = FALSE)
      object <- FindNeighbors(object=object, dims=1:50, reduction = "pca")
      
      ### save PCA coordinates for reproducibility
      pcaCoordinates = object@reductions[['pca']]@cell.embeddings
      cells_PCA = rownames(pcaCoordinates) #CellIds with Sample number as prefix
      pcaCoordinates = as_tibble(pcaCoordinates) #UMAP coordinates in tibble
      pcaCoordinates = pcaCoordinates %>% mutate(cellID = cells_PCA)
      write_csv(pcaCoordinates, file= glue("{pcaCoordinatesDirectory}/{dataset}__{condition}__pcaCoordinates__{dbl_act}.csv"))
      
      # cluster at different resolutions and extract doublet-content information from each cluster
      resolutions =  c(0.2, 0.4, 0.6, 0.8)
      for (resolution in resolutions) {
        object <- FindClusters(object=object, resolution=resolution, algorithm = 1, verbose = FALSE) #clusters based on Louvain algorithm (algorithm = 1)
        cluster_column_name <- glue("RNA_snn_res.{resolution}")
        num_clusters <- length(unique(object@meta.data[[cluster_column_name]]))
        object <- RunUMAP(object = object, reduction = "pca", dims = 1:50, verbose = FALSE)
        
        for (cluster in unique(object@meta.data[[cluster_column_name]])) {
          cluster_cells <- object@meta.data[object@meta.data[[cluster_column_name]] == cluster, ]
          num_doublets <- sum(cluster_cells$label == "doublet", na.rm = TRUE)
          num_singlets <- sum(cluster_cells$label == "singlet", na.rm = TRUE)
          pct_doublets <- mean(cluster_cells$label == "doublet", na.rm = TRUE) * 100  # Assuming 'label' column exists
          
          print(paste("Dataset:", dataset))
          print(paste("Actual Doublets:", dbl_act))
          print(paste("Resolution:", resolution))
          print(paste("Number of Clusters:", num_clusters))
          # print(paste("Cluster:", cluster))
          # print(paste("Number of Doublets:", num_doublets))
          # print(paste("Number of Singlets:", num_singlets))
          print(paste("Percentage of Doublets:", pct_doublets))
          
          # Append the information to the results dataframe
          results <- rbind(results, data.frame(dataset = dataset, 
                                               dbl_act = dbl_act, 
                                               resolution = resolution, 
                                               numClusters = num_clusters, 
                                               cluster = cluster, 
                                               num_doublets = num_doublets, 
                                               num_singlets = num_singlets, 
                                               pct_doublets = pct_doublets))
          
          #print(results)
        }
        
        saveRDS(object, glue("{clusteredObjectsDirectory}/{dataset}_clusteringStabilityObject.rds"))
        
      }
      
    } else {
      message(glue("File {file_path} does not exist."))
    }
  }
  write.csv(results, file = glue("{resultsDirectory}/{dataset}__clusteringStability.csv"), row.names = FALSE)
  print(results)
}




#############

object_joined = JoinLayers(object)
object_joined_norm = NormalizeData(object_joined)
object_joined_norm_var = FindVariableFeatures(object_joined_norm)
object_joined_norm_var_scaled = ScaleData(object_joined_norm_var)
object_joined_sct = SCTransform(object_joined)

object.subset = ScaleData(object.subset)
object.subset = FindVariableFeatures(object.subset)
biorxiv = readRDS(glue("{clusteredObjectsDirectory}/{dataset}_clusteringStabilityObject.rds"))










##########################################################################################################################
# dbscan
##########################################################################################################################


set.seed(23)
#datasets <- list.dirs(functionalObjectsDirectory, full.names = FALSE, recursive = FALSE)

datasets_initial = c("FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "non_cancer","smartseq3_umis", "watermelon")
datasets_overflow  = c("cellTag", "ClonMapper","LARRY", "smartseq3_reads", "SPLINTR", "TREX")

for (dataset in datasets_initial) {
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
      object = JoinLayers(object) # there are many samples which all need to be clustered together. 
      
      ### normalization
      object = NormalizeData(object)
      object = FindVariableFeatures(object)
      object = ScaleData(object)
      
      ### dimension reduction
      object <- RunPCA(object = object, verbose = FALSE)
      object <- FindNeighbors(object=object, dims=1:50, reduction = "pca")
      
      ### save PCA coordinates for reproducibility
      pcaCoordinates = object@reductions[['pca']]@cell.embeddings
      cells_PCA = rownames(pcaCoordinates) #CellIds with Sample number as prefix
      pcaCoordinates = as_tibble(pcaCoordinates) #UMAP coordinates in tibble
      pcaCoordinates = pcaCoordinates %>% mutate(cellID = cells_PCA)
      write_csv(pcaCoordinates, file= glue("{pcaCoordinatesDirectory}/{dataset}__{condition}__pcaCoordinates__{dbl_act}.csv"))
      object <- RunTSNE(object, perplexity=30)
      x <- object@reductions[["tsne"]]@cell.embeddings; dim(x); class(x)
      
      # cluster at different resolutions and extract doublet-content information from each cluster
      resolutions =  c(0.2, 0.4, 0.6, 0.8)
      for (resolution in resolutions) {
        
        res <- dbscan(x, eps = 3, minPts = 5)
        print(table(res$cluster))
        cluster <- names(table(res$cluster)); cluster


        cluster_column_name <- glue("RNA_snn_res.{resolution}")
        num_clusters <- length(unique(res$cluster[res$cluster > 0]))
        #object <- RunUMAP(object = object, reduction = "pca", dims = 1:50, verbose = FALSE)
        
        for (cluster in unique(object@meta.data[[cluster_column_name]])) {
          cluster_cells <- object@meta.data[object@meta.data[[cluster_column_name]] == cluster, ]
          num_doublets <- sum(cluster_cells$label == "doublet", na.rm = TRUE)
          num_singlets <- sum(cluster_cells$label == "singlet", na.rm = TRUE)
          pct_doublets <- mean(cluster_cells$label == "doublet", na.rm = TRUE) * 100  # Assuming 'label' column exists
          
          print(paste("Dataset:", dataset))
          print(paste("Actual Doublets:", dbl_act))
          print(paste("Resolution:", resolution))
          print(paste("Number of Clusters:", num_clusters))
          # print(paste("Cluster:", cluster))
          # print(paste("Number of Doublets:", num_doublets))
          # print(paste("Number of Singlets:", num_singlets))
          print(paste("Percentage of Doublets:", pct_doublets))
          
          # Append the information to the results dataframe
          results <- rbind(results, data.frame(dataset = dataset, 
                                               dbl_act = dbl_act, 
                                               resolution = resolution, 
                                               numClusters = num_clusters, 
                                               cluster = cluster, 
                                               num_doublets = num_doublets, 
                                               num_singlets = num_singlets, 
                                               pct_doublets = pct_doublets))
          
          #print(results)
        }
        
        saveRDS(object, glue("{clusteredObjectsDirectory}/{dataset}_clusteringStabilityObject.rds"))
        
              }
              
            } else {
              message(glue("File {file_path} does not exist."))
            }
          }
          write.csv(results, file = glue("{resultsDirectory}/{dataset}__clusteringStability.csv"), row.names = FALSE)
          print(results)
        }


biorxiv = readRDS(glue("{clusteredObjectsDirectory}/Biorxiv_clusteringStabilityObject.rds"))



