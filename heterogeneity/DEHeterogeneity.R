### assessing heterogeneity using differential expression for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20240213
### Last edited by Madeline E Melzer on 20240222

library(Seurat)
library(glue)
library(purrr)
library(data.table)
library(dplyr)
library(stringr)
library(MAST)
library(devtools)

set.seed(23)

DEDataDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/DE/data"
singletObjectsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/subsampledSingletObjects"


# Function to do DE on a cluster
performDEAnalysis <- function(seurat_obj, cluster, resolution) {
  # Assuming 'seurat_clusters' is the prefix for cluster IDs at different resolutions
  cluster_col_name <- paste0("clusters_res_", as.character(resolution))
  
  # Perform DE analysis using Wilcox, comparing one cluster against all others
  de_results <- FindMarkers(seurat_obj, ident.1 = cluster, group.by = cluster_col_name, test.use = 'wilcox')
  totalGenes = nrow(de_results)
  
  # Filter for significant genes
  significant_genes <- de_results[de_results$p_val_adj < 0.05, ]
  
  # Calculate pUp, pDown, fU, and fL
  pUp <- sum(significant_genes$avg_log2FC > 0) / totalGenes
  pDown <- sum(significant_genes$avg_log2FC < 0) / totalGenes
  fL <- min(significant_genes$avg_log2FC)  # Minimum log fold change for fL
  fU <- max(significant_genes$avg_log2FC)  # Maximum log fold change for fU
  
  return(data.frame(Cluster = cluster, Resolution = resolution, pUp = pUp, pDown = pDown, fU = fU, fL = fL))
}




# Load each Seurat object (sample)
singletFiles <- list.files(singletObjectsDirectory, pattern = "\\.rds$", full.names = TRUE)
sample_names <- gsub(pattern = ".*/|\\.rds$", replacement = "", x = sample_files)

# Initialize an empty data frame to store comprehensive results
results <- data.frame()

for (filepath in singletFiles) {
  print(filepath) # to monitor progress
  parts <- strsplit(basename(filepath), "__")[[1]]
  dataset <- parts[1]
  sample <- gsub("_subsampledSinglets\\.rds$", "", parts[2])
  
  # Load the Seurat object
  seuratObject <- readRDS(filepath)
  
  for (i in c(0.2, 0.4, 0.6, 0.8)) {
    resolution = i
    # Get unique clusters at this resolution
    cluster_col_name <- paste0("clusters_res_", as.character(resolution))
    clusters <- unique(seuratObject@meta.data[[cluster_col_name]])
    
    # Loop through each cluster
    for (cluster in clusters) {
      # Perform DE analysis and collect results
      de_result <- performDEAnalysis(seuratObject, cluster, resolution)
      
      # Add dataset, sample, and cluster information to the results
      de_result$Dataset <- dataset
      de_result$Sample <- sample
      de_result$Cluster <- cluster
      de_result$Resolution <- resolution
      
      # Append to the comprehensive results data frame
      #print(de_result)
      results <- rbind(results, de_result)
      #print(results)
    }
  }
}

# Save the comprehensive DE results to a file
write_csv(results, glue("{DEDataDir}/DEResults_wilcox.csv"))


### troubleshooting
# 
# filepath = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/subsampledSingletObjects/cellTag__B4D21-RNA-r2-1__subsampledSinglets.rds"    
# seuratObject <- readRDS(filepath)
# resolution = 0.6
# # Get unique clusters at this resolution
# cluster_col_name <- paste0("clusters_res_", as.character(resolution))
# clusters <- unique(seuratObject@meta.data[[cluster_col_name]])
# cluster = 7
# 
# # de
# seurat_obj = seuratObject
# # Perform DE analysis using MAST, comparing one cluster against all others
# de_results <- FindMarkers(seurat_obj, ident.1 = cluster, group.by = cluster_col_name, test.use = 'wilcox')
# 
# # Filter for significant genes
# significant_genes <- de_results[de_results$p_val_adj < 0.05, ]
# 
# # Calculate pUp, pDown, fU, and fL
# pUp <- mean(significant_genes$avg_log2FC > 0)
# pDown <- mean(significant_genes$avg_log2FC < 0)
# fL <- min(significant_genes$avg_log2FC)  # Minimum log fold change for fL
# fU <- max(significant_genes$avg_log2FC)  # Maximum log fold change for fU


















s_root <- "/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset"
# output_root <- "/projects/p31666/zzhang/doublet-bchmk/data/DE_output_supp_run"
output_root <- "/projects/p31666/zzhang/doublet-bchmk/data/DE_output"
# TODO: Change if rerunning all
# all_dataset_dir <- list.dirs(s_root, recursive = F)
all_dataset_dir <- c(
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/SPLINTR",
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/TREX"
)


for(cur_dataset_dir in all_dataset_dir){
  cur_dataset_id <- basename(cur_dataset_dir)
  if(file.exists(glue("{cur_dataset_dir}/data/forty.rds")) == FALSE){
    next
  }
  forty_seu <- readRDS(glue("{cur_dataset_dir}/data/forty.rds"))
  twenty_seu <- readRDS(glue("{cur_dataset_dir}/data/twenty.rds"))
  ten_seu <- readRDS(glue("{cur_dataset_dir}/data/ten.rds"))
  control_seu <- readRDS(glue("{cur_dataset_dir}/data/control.rds"))
  cur_dataset_all_s <- list(
    "forty" = forty_seu,
    "twenty" = twenty_seu,
    "ten" = ten_seu,
    "control" = control_seu
  )
  Idents(forty_seu) <- "seurat_clusters"
  cur_clusters <- Idents(forty_seu)|>levels()
  combinations <- combn(cur_clusters, 2)
  num_combinations_to_sample <- min(3, ncol(combinations))
  selected_combinations <- combinations[, sample(ncol(combinations), num_combinations_to_sample)]

  cur_dataset_out_dir <- glue("{output_root}/{seed}/{cur_dataset_id}")
  dir.create(cur_dataset_out_dir, recursive = TRUE, showWarnings = FALSE)
  for(cur_dataset_dbl_pct in names(cur_dataset_all_s)){
    cur_s <- cur_dataset_all_s[[cur_dataset_dbl_pct]]
    cells_to_keep <- sample(Cells(cur_s), 100)
    cur_s <- subset(cur_s, cells = cells_to_keep)
    
    cur_s[["RNA"]] <- JoinLayers(cur_s[["RNA"]])
    Idents(cur_s) <- "seurat_clusters"
    for(cur_col in ncol(selected_combinations)){
      cur_clus_1 <- selected_combinations[1, cur_col]
      cur_clus_2 <- selected_combinations[2, cur_col]

      cur_out_file <- glue("{cur_dataset_out_dir}/{cur_dataset_dbl_pct}_{cur_clus_1}_vs._{cur_clus_2}__wilcox.csv")
      if(file.exists(cur_out_file) == FALSE){
        cur_DE <- FindMarkers(cur_s, ident.1 = cur_clus_1, ident.2 = cur_clus_2)
        cur_DE[["gene"]] <- rownames(cur_DE)
        write.csv(cur_DE, file = cur_out_file, quote = F, row.names = F)
      }else{
        print(glue("Skip: {cur_out_file}"))
      }

      cur_out_file <- glue("{cur_dataset_out_dir}/{cur_dataset_dbl_pct}_{cur_clus_1}_vs._{cur_clus_2}__MAST.csv")
      if(file.exists(cur_out_file) == FALSE){
        cur_DE <- FindMarkers(cur_s, ident.1 = cur_clus_1, ident.2 = cur_clus_2, test.use = "MAST")
        cur_DE[["gene"]] <- rownames(cur_DE)
        write.csv(cur_DE, file = cur_out_file, quote = F, row.names = F)
      }else{
        print(glue("Skip: {cur_out_file}"))
      }
    }
  }
}