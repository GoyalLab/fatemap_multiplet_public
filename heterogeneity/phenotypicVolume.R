### assessing phenotypic-volume based heterogeneity for Zhang Melzer et al 2024
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
