#Adopting TREX (RatzEtAl Nat Neuro 2022) data for ZhangMelzerEtAl Cell Genom 2023
#Last updated 20231023 by Madeline E Melzer

library(tidyverse)
library(DropletUtils)
library(Seurat)
library(ggplot2)


dataDirectory = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/"

brain1_39.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644058_10x39_filtered_gene_bc_matrix/")
brain1_40.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644059_10x40_filtered_gene_bc_matrix/")
brain1_41.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644060_10x41_filtered_gene_bc_matrix/")


brain1_39 = CreateSeuratObject(brain1_39.data)
brain1_40 = CreateSeuratObject(brain1_40.data)
brain1_41 = CreateSeuratObject(brain1_40.data)

brain1 = merge(brain1_39, y = c(brain1_40, brain1_41), add.cell.ids = c("10x39", "10x40", "10x41"))
#saveRDS(brain1, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain1.rds")
brain1_counts_matrix <- as(brain1@assays$RNA@counts, "sparseMatrix") # Ensure the assay name is correct
DropletUtils::write10xCounts(
  path = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain1_10X_01/",
  x = brain1_counts_matrix
)

brain2_42.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644061_10x42_filtered_gene_bc_matrix/")
brain2_43.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644062_10x43_filtered_gene_bc_matrix/")
brain2_44.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644063_10x44_filtered_gene_bc_matrix/")

brain2_42 = CreateSeuratObject(brain2_42.data)
brain2_43 = CreateSeuratObject(brain2_43.data)
brain2_44 = CreateSeuratObject(brain2_44.data)

brain2 = merge(brain2_42, y = c(brain2_43, brain2_44), add.cell.ids = c("10x42", "10x43", "10x44"))
#saveRDS(brain2, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain2.rds")
brain2_counts_matrix <- as(brain2@assays$RNA@counts, "sparseMatrix") # Ensure the assay name is correct
DropletUtils::write10xCounts(
  path = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain2_10X/",
  x = brain2_counts_matrix
)


brain3_52.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644064_10x52_filtered_gene_bc_matrix/")
brain3_53.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644065_10x53_filtered_gene_bc_matrix/")
brain3_54.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644066_10x54_filtered_gene_bc_matrix/")

brain3_52 = CreateSeuratObject(brain3_52.data)
brain3_53 = CreateSeuratObject(brain3_53.data)
brain3_54 = CreateSeuratObject(brain3_54.data)

brain3 = merge(brain3_52, y = c(brain3_53, brain3_54), add.cell.ids = c("10x52", "10x53", "10x54"))
#saveRDS(brain3, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain3.rds")
brain3_counts_matrix <- as(brain3@assays$RNA@counts, "sparseMatrix") # Ensure the assay name is correct
DropletUtils::write10xCounts(
  path = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain3_10X/",
  x = brain3_counts_matrix
)

brain4_65.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644067_10x65_filtered_gene_bc_matrix/")
brain4_66.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644068_10x66_filtered_gene_bc_matrix/")
brain4_67.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/GSM4644069_10x67_filtered_gene_bc_matrix/")

brain4_65 = CreateSeuratObject(brain4_65.data)
brain4_66 = CreateSeuratObject(brain4_66.data)
brain4_67 = CreateSeuratObject(brain4_67.data)

brain4 = merge(brain4_65, y = c(brain4_66, brain4_67), add.cell.ids = c("10x65", "10x66", "10x67"))
#saveRDS(brain4, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain4.rds")
brain4_counts_matrix <- as(brain4@assays$RNA@counts, "sparseMatrix") # Ensure the assay name is correct
DropletUtils::write10xCounts(
  path = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/TREX/GSE153424_RAW/brain4_10X/",
  x = brain4_counts_matrix
)








