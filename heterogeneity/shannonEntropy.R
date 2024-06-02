### using shannon entropy for assessing heterogeneity for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20231220
### Last edited by Madeline E Melzer on 20240222

library(tidyverse)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(glue)
library(umap)
library(vegan)
library(ggpubr)
library(permute)
options(future.globals.maxSize = 4000 * 1024^2)

set.seed = 23

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
# calculating Shannon Entropy
##########################################################################################################################
# adapted from Goyal et al. 2023, FM01_ShannonEntropy_barcodes.R  and Fennell, Vassiliadis et al., 2022, scRNA_transcriptional_diversity_analysis.Rmd

shannonEntropy <- data.frame()
singletObjectsDirectory = "/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/singletObjects/subsampledSingletObjects"
singletFiles <- list.files(singletObjectsDirectory, pattern = "\\.rds$", full.names = TRUE)

for(filepath in singletFiles) {
  # Extract dataset and sample from the file name
  print(filepath)
  parts <- strsplit(basename(filepath), "__")[[1]]
  dataset <- parts[1]
  sample <- gsub("_subsampledSinglets\\.rds$", "", parts[2])
  
  seuratObject <- readRDS(filepath)
  
  for (i in c(0.2, 0.4, 0.6, 0.8)) {
    resolution = i
    seurat_clusters <- paste0("clusters_res_", resolution)
    cluster_counts <- table(seuratObject@meta.data[[seurat_clusters]])
    cluster_matrix <- matrix(cluster_counts, ncol = 1)
    shannon_diversity_cluster <- diversity(cluster_matrix, index = "shannon")
    numClusters <- length(unique(seuratObject@meta.data[[seurat_clusters]]))
    
    ### Shuffling cluster numbers as a control
    cluster_data <- seuratObject[[paste0("clusters_res_", as.character(resolution))]]
    cluster_data_shuff = cluster_data
    cluster_data_shuff[[1]] <- sample(cluster_data[[1]])
    seuratObject@meta.data[, paste0("clusters_res_", as.character(resolution), "_shuffled")] <- cluster_data_shuff
    
    shuffled_clusters <- paste0("clusters_res_", as.character(resolution), "_shuffled")
    shuffled_counts <- table(seuratObject@meta.data[[shuffled_clusters]])
    shuffled_matrix <- matrix(shuffled_counts, ncol = 1)
    shannon_diversity_shuffled <- diversity(shuffled_matrix, index = "shannon")
    
    ### uniform cluster distribution as a control (maximum equitability)
    uniform_counts <- rep(1, numClusters)
    shannon_diversity_uniform <- diversity(matrix(uniform_counts, ncol = 1), index = "shannon")
    
    temp_df <- data.frame(
      dataset = dataset,
      sample = sample,
      resolution = resolution,
      numClusters = numClusters,
      clusterDiversity = shannon_diversity_cluster,
      shuffledDiversity = shannon_diversity_shuffled,
      uniformDiversity = shannon_diversity_uniform
    )
    shannonEntropy <- rbind(shannonEntropy, temp_df)
  }
}

write_csv(shannonEntropy, "~/ZhangMelzerEtAl/data/shannonEntropy/shannonEntropy_2.csv")


