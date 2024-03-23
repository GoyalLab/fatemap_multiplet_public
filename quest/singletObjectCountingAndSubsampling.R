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
    if(numSinglets > 1000) {
      subsampledObject <- getSubsampledClusteredSeurat(seuratObject)  # Assuming getSubsampledSeurat is defined elsewhere
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
  analysis.seurat = seuratObject
  cells <- sample(colnames(analysis.seurat), size =numCells, replace=F)
  analysis.seurat.subset <- analysis.seurat[, cells]
  
  ### normalization
  analysis.seurat.subset = NormalizeData(analysis.seurat.subset)
  analysis.seurat.subset = FindVariableFeatures(analysis.seurat.subset)
  analysis.seurat.subset = ScaleData(analysis.seurat.subset)
  
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



##########################################################################################################################
# merge Jain et al. ("Biorxiv") DMSO A and B to make it have over 1000 cells (only 890 cells in smallest sample)
##########################################################################################################################

sample1_dir = "/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/Biorxiv/10X/1_DMSO_A"
sample2_dir = "/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/Biorxiv/10X/2_DMSO_B"

dataset = "Biorxiv"
sample = "1_DMSO_A"
#setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/", dataset, "/fatemapID/")
#sample1 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.01)
#sample1 = Read10X(data_dir)
#sample1 = CreateSeuratObject(sample1)
sample1@meta.data$sample <- "1_DMSO_A"


sample = "2_DMSO_B"
#setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/fatemapID/")
#sample2 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.01)
sample2 = Read10X(data_dir)
sample2 = CreateSeuratObject(sample2)
sample2@meta.data$sample <- "2_DMSO_B"

#saveRDS(sample1, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/sample1/sample1.rds")
#saveRDS(sample2, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/sample2/sample2.rds")

sample1 = readRDS("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/Biorxiv/1_DMSO_A_singletsOnly.rds")
sample1@meta.data$sample <- "1_DMSO_A"

sample2 = readRDS("/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/Biorxiv/2_DMSO_B_singletsOnly.rds")
sample2@meta.data$sample <- "2_DMSO_B"

sample1 = PercentageFeatureSet(sample1, pattern = "^MT-", col.name = "percent.mt")
sample2 = PercentageFeatureSet(sample2, pattern = "^MT-", col.name = "percent.mt")

s1s2_scTransform <- merge(sample1, y = sample2, project = "S1S2")
s1s2_scTransform.joined <- JoinLayers(s1s2_scTransform) #this is because of the 2 layers/
s1s2_scTransform.list <- SplitObject(s1s2_scTransform.joined, split.by = "sample")
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



###### getting clustered, subsampled seurat from this data instead of re-normalizing, etc.
dataset = "Biorxiv"
numCells = 1000
#Biorxiv_subsampled = getSubsampledClusteredSeurat(s1s2_scTransform.integrated, 1000)

### subsample numCells cells from the sample
analysis.seurat = s1s2_scTransform.integrated
cells <- sample(colnames(analysis.seurat), size =numCells, replace=F)
analysis.seurat.subset <- analysis.seurat[, cells]

### normalization
# analysis.seurat.subset = NormalizeData(analysis.seurat.subset)
# analysis.seurat.subset = FindVariableFeatures(analysis.seurat.subset)
# analysis.seurat.subset = ScaleData(analysis.seurat.subset)

### dimension reduction
# analysis.seurat.subset <- RunPCA(object = analysis.seurat.subset, verbose = FALSE)
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

saveRDS(analysis.seurat.subset, file = glue("{subsampledSingletObjectsDirectory}/{dataset}__{sample}__subsampledSinglets.rds"))





