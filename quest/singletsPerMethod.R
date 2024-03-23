### exteacting the singlets identified by each doublet detection method for Zhang, Melzer et al 2024
### Created by Madeline E Melzer on 20240211
### Last edited by Madeline E Melzer on 20240316
# last edit was to get all of the act_0.08 data so that Charles can use it for TNR

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
options(future.globals.maxSize = 4000 * 1024^2)

doubletObjectsDirectory = "/projects/b1042/GoyalLab/zzhang/doublet_objects"
doubletInfoDirectory = "~/ZhangMelzerEtAl/data/doubletInfo0.08act"

##########################################################################################################################
# doublet_finder
##########################################################################################################################

detectionMethod = "doublet_finder"
dblMethodPath <- file.path(doubletObjectsDirectory, detectionMethod)
# Get the list of dbl_act directories
dblActs <- list.dirs(dblMethodPath, full.names = FALSE, recursive = FALSE)

# Iterate over dbl_act
for (dbl_act in dblActs) {
  dblActPath <- file.path(dblMethodPath, dbl_act)
  
  # Get the list of .rds files
  rdsFiles <- list.files(dblActPath, pattern = "\\.rds$", full.names = TRUE)
  
  # Iterate over .rds files
  for (file in rdsFiles) {
    filename <- basename(file)
    filename <- tools::file_path_sans_ext(filename)
    if (grepl("exp_([0-9]+(?:\\.[0-9]+)?)__act_0\\.08\\b", filename)) {
      seuratObject = readRDS(file)
      
      metadata <- seuratObject@meta.data %>%
        rownames_to_column("cellID") %>%
        select(cellID, nCount_RNA, nFeature_RNA, label, contains("pAnn"), contains("DF.classifications"))
      
      write_csv(metadata, file = glue("~/ZhangMelzerEtAl/data/doubletInfo_0.08act/{detectionMethod}/{filename}___doubletData.csv"))
    }
  }
}

##########################################################################################################################
# doublet_cell
##########################################################################################################################

detectionMethod = "doublet_cell"
dblMethodPath <- file.path(doubletObjectsDirectory, detectionMethod)
# Get the list of dbl_act directories
dblActs <- list.dirs(dblMethodPath, full.names = FALSE, recursive = FALSE)

# Iterate over dbl_act
for (dbl_act in dblActs) {
  dblActPath <- file.path(dblMethodPath, dbl_act)
  
  # Get the list of .rds files
  rdsFiles <- list.files(dblActPath, pattern = "\\.rds$", full.names = TRUE)
  
  # Iterate over .rds files
  for (file in rdsFiles) {
    filename <- basename(file)
    filename <- tools::file_path_sans_ext(filename)
    
    if (grepl("exp_([0-9]+(?:\\.[0-9]+)?)__act_0\\.08\\b", filename)) {
      sce = readRDS(file)
      
      metadata <- as.data.frame(colData(sce)) %>%
        rownames_to_column("cellID") %>%
        select(cellID, nCount_RNA, nFeature_RNA, label, scDblFinder.class, scDblFinder.score, scDblFinder.weighted, scDblFinder.cxds_score)
      
      write_csv(metadata, file = glue("~/ZhangMelzerEtAl/data/doubletInfo_0.08act/{detectionMethod}/{filename}___doubletData.csv"))
    }
  }
}


##########################################################################################################################
# hybrid
##########################################################################################################################
detectionMethod = "hybrid"
dblMethodPath <- file.path(doubletObjectsDirectory, detectionMethod)
# Get the list of dbl_act directories
dblActs <- list.dirs(dblMethodPath, full.names = FALSE, recursive = FALSE)

# Iterate over dbl_act
for (dbl_act in dblActs) {
  dblActPath <- file.path(dblMethodPath, dbl_act)
  
  # Get the list of .rds files
  rdsFiles <- list.files(dblActPath, pattern = "\\.rds$", full.names = TRUE)
  
  # Iterate over .rds files
  for (file in rdsFiles) {
    filename <- basename(file)
    filename <- tools::file_path_sans_ext(filename)
    
    if (grepl("exp_([0-9]+(?:\\.[0-9]+)?)__act_0\\.08\\b", filename)) {
      sce = readRDS(file)
      
      metadata <- as.data.frame(colData(sce)) %>%
        rownames_to_column("cellID") %>%
        select(cellID, nCount_RNA, nFeature_RNA, label, hybrid_score, hybrid_call, cxds_score, cxds_call, bcds_score, bcds_call)
      
      write_csv(metadata, file = glue("~/ZhangMelzerEtAl/data/doubletInfo_0.08act/{detectionMethod}/{filename}___doubletData.csv"))
    }
  }
}







##########################################################################################################################
# scrublet (do not need because already .csv files)
##########################################################################################################################

rds = readRDS("/projects/b1042/GoyalLab/zzhang/doublet_objects/doublet_cell/act__0.1/Biorxiv___1_DMSO_A___exp_0.1__act_0.1.rds")
nrow(colData(rds))

file = "/projects/b1042/GoyalLab/zzhang/doublet_objects/doublet_finder/act__0.1/Biorxiv___1_DMSO_A___exp_0.1__act_0.1.rds"

### testing i think 
doubletFinder = readRDS("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doublet_objects/doublet_finder/act__0.1/Biorxiv___1_DMSO_A___exp_0.1__act_0.1.rds")
seuratObject <- doubletFinder

detectionMethod = "doublet_cell"
file = glue("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doublet_objects/{detectionMethod}/act__0.1/Biorxiv___1_DMSO_A___exp_0.1__act_0.1.rds")
seuratObject = readRDS(file)
filename <- basename(file)
filename <- tools::file_path_sans_ext(filename)

metadata <- seuratObject@colData %>%
  rownames_to_column("cellID") %>%
  select(cellID, nCount_RNA, nFeature_RNA, label, scDblFinder.class, scDblFinder.score, )

write_csv(metadata, file = glue("{doubletInfoDirectory}/{detectionMethod}/{filename}___doubletData.csv"))
