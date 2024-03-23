### counting and subsampling 1000 singlets for heterogeneity analysis for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20240222
### Last edited by Madeline E Melzer on 20240222

library(tidyverse)
library(dplyr)
library(Seurat)
library(glue)
library(readr)
library(umap)

singletObjectsDirectory = "/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/"
subsampledSingletObjectsDirectory = "/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/subsampledSingletObjects"

##########################################################################################################################
# counting singlets per object, subsampling 1000 cells (if >1000 cells), and saving object in subsampledSingletObjects directory and singlet counts in singletCountBySample.csv
##########################################################################################################################

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

singletCounts <- data.frame(dataset = character(), sample = character(), numSinglets = integer())

# Loop through each dataset and process the Seurat singlet-only objects
for(dataset in datasets) {
  datasetDir <- file.path(singletObjectsDirectory, dataset)
  seuratFiles <- list.files(datasetDir, pattern = "_singletsOnly\\.rds$", full.names = TRUE)
  
  for(filepath in seuratFiles) {
    
    sample <- gsub("(.*)_singletsOnly\\.rds$", "\\1", basename(filepath)) # Extract the sample name from the file name
    seuratObject <- readRDS(filepath) # Load the Seurat object
    numSinglets <- ncol(seuratObject) # Count the cells
    singletCounts <- rbind(singletCounts, data.frame(dataset = dataset, sample = sample, numSinglets = numSinglets)) # Add to the data frame
    
    # Check if cell count is more than 1000 and subsample if it is
    if(cellCount > 1000) {
      subsampledObject <- getSubsampledSeurat(seuratObject)  # Assuming getSubsampledSeurat is defined elsewhere
      saveRDS(subsampledObject, file = glue("{subsampledSingletObjectsDirectory}/{dataset}__{sample}__subsampledSinglets.rds"))
    }
  }
}

# Save the singlet counts to a CSV file
write_csv(singletCounts, glue("{singletObjectsDirectory}/singletCountBySample.csv"))




##########################################################################################################################
# FUNCTIONS
##########################################################################################################################

getSubsampledClusteredSeurat <- function (seuratObject, numCells) {
  
  ### subsample numCells cells from the sample
  analysis.seurat = object
  cells <- sample(colnames(analysis.seurat), size =numCells, replace=F)
  analysis.seurat.subset <- analysis.seurat[, cells]
  
  ### normalization
  object = NormalizeData(object)
  object = FindVariableFeatures(object)
  object = ScaleData(object)
  
  ### dimension reduction
  analysis.seurat.subset <- RunPCA(object = analysis.seurat.subset, verbose = FALSE)
  analysis.seurat.subset <- FindNeighbors(object=analysis.seurat.subset, dims=1:50, reduction = "pca", force.recalc = TRUE)
  
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

