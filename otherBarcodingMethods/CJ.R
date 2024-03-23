#Adopting data from CJE (Jeff Mold, biorXiv, 2022) data for ZhangMelzerEtAl 2024
#Created by Madeline E Melzer on 20231030
# Last updated 20240207 by Madeline E Melzer

library(tidyverse)
library(DropletUtils)
library(Seurat)
library(ggplot2)
library(dplyr)

#this comes from the Smartseq3.dgecounts.rds provided by CJE, extracted using "Smartseq3.dgecounts$umicount$inex$all" and saved with DropletUtils::write10xCounts()
data = readRDS("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/smartseq3umiinexallcounts.rds")
data <- data[rownames(data) != "EGFP-30N", ] #do not want this in final counts matrix

lineageBC_umis = read_tsv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/providedByCJ/lineageBC_results_UMIonly.txt")
lineageBC_reads = read_tsv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/providedByCJ/lineageBC_results.txt")

brain1_readCountMatrix = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/providedByCJ/Archive/brain1/read_count_matrix.csv")
brain2_readCountMatrix = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/providedByCJ/Archive/brain2/read_count_matrix.csv")

brain1_cells = colnames(brain1_readCountMatrix)[2:ncol(brain1_readCountMatrix)] #the first column is the title of the barcodes column, "...1", and we do not want this
brain2_cells = colnames(brain2_readCountMatrix)[2:ncol(brain2_readCountMatrix)] #the first column is the title of the barcodes column, "...1", and we do not want this

#get umi counts for each cell-barcode combination in each brain
brain1_umis <- lineageBC_umis %>%
  dplyr::filter(BC %in% brain1_cells)
brain2_umis <- lineageBC_umis %>%
  dplyr::filter(BC %in% brain2_cells)

#expand rows based on umi count
brain1_umis_expanded <- brain1_umis %>%
  uncount(weights = num_umis, .remove = TRUE) %>%
  mutate(sample = "brain1")
brain2_umis_expanded <- brain2_umis %>%
  uncount(weights = num_umis, .remove = TRUE) %>%
  mutate(sample = "brain2")

#rename, combine, and save for umi counts
umis <- bind_rows(brain1_umis_expanded, brain2_umis_expanded)
umis = umis %>% dplyr::rename(barcode = lineageBC,
                       cellID = BC)
write.csv(umis, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/smartseq3_umis/all_brains_cellID_barcode_UMI.csv", row.names = FALSE)
length(unique(umis$cellID))

#get read counts for each cell-barcode combination in each brain
brain1_reads <- lineageBC_reads %>%
  dplyr::filter(BC %in% brain1_cells)
brain2_reads <- lineageBC_reads %>%
  dplyr::filter(BC %in% brain2_cells)

#expand rows based on read count
brain1_reads_expanded <- brain1_reads %>%
  uncount(weights = num_reads, .remove = TRUE) %>%
  mutate(sample = "brain1")

brain2_reads_expanded <- brain2_reads %>%
  uncount(weights = num_reads, .remove = TRUE) %>%
  mutate(sample = "brain2")

#rename, combine, and save for read counts
reads <- bind_rows(brain1_reads_expanded, brain2_reads_expanded)
reads = reads %>% dplyr::rename(barcode = lineageBC,
                              cellID = BC)
write.csv(reads, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/smartseq3_reads/all_brains_cellID_barcode_read.csv", row.names = FALSE)


# Subset the original counts matrix to only include the cells for a brain
brain1_data <- data[, brain1_cells]
print(length(brain1_cells[brain1_cells %in% colnames(brain1_data)])) #checking overlap, 1802
DropletUtils::write10xCounts(
  path = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/10X/brain1/",
  x = brain1_data
)

brain2_data <- data[, brain2_cells]
print(length(brain2_cells[brain2_cells %in% colnames(brain2_data)])) #checking overlap, 992
DropletUtils::write10xCounts(
  path = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/10X/brain2/",
  x = brain2_data
)

# formatting check: making sure I can load it into Seurat

brain1_10X = Read10X(data.dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/10X/brain1/")
brain2_10X = Read10X(data.dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/CJ/smartseq3/10X/brain2/")
#note- renamed to features.tsv for consistency with other method names

brain1_seurat = CreateSeuratObject(counts = brain1_10X, min.cells = 3, min.features = 200)
colnames(brain1_seurat@assays$RNA)
rownames(brain1_seurat@assays$RNA)

brain2_seurat = CreateSeuratObject(counts = brain2_10X, min.cells = 3, min.features = 200)
colnames(brain2_seurat@assays$RNA)
rownames(brain2_seurat@assays$RNA)
