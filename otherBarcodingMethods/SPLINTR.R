#Adopting SPLINTR (FennelandVassiliadisEtAl Nature 2022) data for ZhangMelzerEtAl 2023
#Created by Madeline E Melzer on 20231025
# Last updated 20231025 by Madeline E Melzer

library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
library(ggplot2)

########################## Demultiplexing hashed, pooled cells for chemotherapy treated datasets ################################################

pool1.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_1/", strip.suffix = TRUE)
pool1HTO.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_1_HTO/", gene.column = 1)

#pool1HTO.matrix = paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_1_HTO/matrix.mtx.gz")
pool1.joint.bcs = intersect(colnames(pool1.data), colnames(pool1HTO.data))

pool1.data = pool1.data[, pool1.joint.bcs, drop = FALSE]
pool1HTO.data = as.matrix(pool1HTO.data[, pool1.joint.bcs, drop = FALSE])

rownames(pool1HTO.data) #check if row names are the HTO names/types

# Setup Seurat object
pool1.hashtag = CreateSeuratObject(counts = pool1.data)
# Normalize RNA data with log normalization
pool1.hashtag <- NormalizeData(pool1.hashtag)
# Find and scale variable features
pool1.hashtag <- FindVariableFeatures(pool1.hashtag, selection.method = "mean.var.plot")
pool1.hashtag <- ScaleData(pool1.hashtag, features = VariableFeatures(pool1.hashtag))

# Add HTO data as a new assay independent from RNA
pool1.hashtag[["HTO"]] = CreateAssayObject(counts = pool1HTO.data)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pool1.hashtag <- NormalizeData(pool1.hashtag, assay = "HTO", normalization.method = "CLR")

pool1.hashtag.tagged = pool1.hashtag

# Assuming `features_to_remove` is a character vector of the row names you want to remove from the "HTO" assay
features_to_remove <- rownames(pool1.hashtag.tagged[["HTO"]])[c(4, 6, 7)]

features_to_remove <- c("TotalSeq-A0304-AAAGCATTCTTCACG", "TotalSeq-A0306-TATGCTGCCACGGTA", "unmapped")

features_to_remove_idx <- which(rownames(pool1.hashtag.tagged[["HTO"]]) %in% features_to_remove)

pool1.hashtag.tagged[["HTO"]]@counts <- pool1.hashtag.tagged[["HTO"]]@counts[-features_to_remove_idx, , drop = FALSE]
pool1.hashtag.tagged[["HTO"]]@data <- pool1.hashtag.tagged[["HTO"]]@data[-features_to_remove_idx, , drop = FALSE]
pool1.hashtag.tagged[["HTO"]]@meta.features <- pool1.hashtag.tagged[["HTO"]]@meta.features[-features_to_remove_idx, , drop = FALSE]
# Remove features
#pool1.hashtag.tagged[["HTO"]]@counts <- pool1.hashtag.tagged[["HTO"]]@counts[-which(rownames(pool1.hashtag.tagged[["HTO"]]) %in% features_to_remove), ]
#pool1.hashtag.tagged[["HTO"]]@data <- pool1.hashtag.tagged[["HTO"]]@data[-which(rownames(pool1.hashtag.tagged[["HTO"]]) %in% features_to_remove), ]
#pool1.hashtag.tagged[["HTO"]] = pool1.hashtag.tagged[["HTO"]][-which(rownames(pool1.hashtag.tagged[["HTO"]]) %in% features_to_remove), ]

#rownames(pool1.hashtag.tagged) <- c("VEH-M15-BM", "VEH-M05-BM", "CHEMO-M10-BM", "CHEMO-M02-BM")

pool1.hashtag.demux <- HTODemux(pool1.hashtag.tagged, assay = "HTO", positive.quantile = 0.99)

table(pool1.hashtag.demux$HTO_classification.global)
Idents(pool1.hashtag.demux) <- "HTO_classification.global"
VlnPlot(pool1.hashtag.demux, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Convert the named vector to a data frame
hto_maxID_df <- data.frame(
  CellName = paste0(names(pool1.hashtag.demux$HTO_maxID), "-1"),
  HTO_maxID = pool1.hashtag.demux$HTO_maxID,
  HTO_classification = pool1.hashtag.demux$HTO_classification,
  HTO_classification.global = pool1.hashtag.demux$HTO_classification.global,
  hash.ID = pool1.hashtag.demux$hash.ID,
  row.names = NULL
)
#write.csv(hto_maxID_df, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_1_HTO_maxID.csv", row.names = FALSE)



####### hashing pool 2 ################


pool2.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_2/", strip.suffix = TRUE)
pool2HTO.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_2_HTO/", gene.column = 1)

#pool2HTO.matrix = paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_1_HTO/matrix.mtx.gz")
pool2.joint.bcs = intersect(colnames(pool2.data), colnames(pool2HTO.data))

pool2.data = pool2.data[, pool2.joint.bcs, drop = FALSE]
pool2HTO.data = as.matrix(pool2HTO.data[, pool2.joint.bcs, drop = FALSE])

rownames(pool2HTO.data) #check if row names are the HTO names/types

# Setup Seurat object
pool2.hashtag = CreateSeuratObject(counts = pool2.data)
# Normalize RNA data with log normalization
pool2.hashtag <- NormalizeData(pool2.hashtag)
# Find and scale variable features
pool2.hashtag <- FindVariableFeatures(pool2.hashtag, selection.method = "mean.var.plot")
pool2.hashtag <- ScaleData(pool2.hashtag, features = VariableFeatures(pool2.hashtag))

# Add HTO data as a new assay independent from RNA
pool2.hashtag[["HTO"]] = CreateAssayObject(counts = pool2HTO.data)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pool2.hashtag <- NormalizeData(pool2.hashtag, assay = "HTO", normalization.method = "CLR")

pool2.hashtag.tagged = pool2.hashtag

# Assuming `features_to_remove` is a character vector of the row names you want to remove from the "HTO" assay
features_to_remove <- rownames(pool2.hashtag.tagged[["HTO"]])[c(4, 6, 7)]

features_to_remove <- c("TotalSeq-A0304-AAAGCATTCTTCACG", "TotalSeq-A0306-TATGCTGCCACGGTA", "unmapped")  # Replace with actual feature names

features_to_remove_idx <- which(rownames(pool2.hashtag.tagged[["HTO"]]) %in% features_to_remove)

pool2.hashtag.tagged[["HTO"]]@counts <- pool2.hashtag.tagged[["HTO"]]@counts[-features_to_remove_idx, , drop = FALSE]
pool2.hashtag.tagged[["HTO"]]@data <- pool2.hashtag.tagged[["HTO"]]@data[-features_to_remove_idx, , drop = FALSE]
pool2.hashtag.tagged[["HTO"]]@meta.features <- pool2.hashtag.tagged[["HTO"]]@meta.features[-features_to_remove_idx, , drop = FALSE]
# Remove features
#pool2.hashtag.tagged[["HTO"]]@counts <- pool2.hashtag.tagged[["HTO"]]@counts[-which(rownames(pool2.hashtag.tagged[["HTO"]]) %in% features_to_remove), ]
#pool2.hashtag.tagged[["HTO"]]@data <- pool2.hashtag.tagged[["HTO"]]@data[-which(rownames(pool2.hashtag.tagged[["HTO"]]) %in% features_to_remove), ]
#pool2.hashtag.tagged[["HTO"]] = pool2.hashtag.tagged[["HTO"]][-which(rownames(pool2.hashtag.tagged[["HTO"]]) %in% features_to_remove), ]

#rownames(pool2.hashtag.tagged) <- c("VEH-M15-BM", "VEH-M05-BM", "CHEMO-M10-BM", "CHEMO-M02-BM")

pool2.hashtag.demux <- HTODemux(pool2.hashtag.tagged, assay = "HTO", positive.quantile = 0.99)

table(pool2.hashtag.demux$HTO_classification.global)
Idents(pool2.hashtag.demux) <- "HTO_classification.global"
VlnPlot(pool2.hashtag.demux, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

############### Labeling hashed cell singlets properly ####################################

POOL_1 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_2_HTO_maxID.csv")
POOL_1_map = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_1_map.csv")

POOL_1 = POOL_1 %>%
  rename(Cell.10X.Barcode = CellName) %>%
  inner_join(POOL_1_map, by = "Cell.10X.Barcode") %>%
  filter(!(hash.ID %in% c("Doublet", "Negative"))) %>%
  mutate(hash.ID = case_when(
    hash.ID == "TotalSeq-A0301-ACCCACCAGTAAGAC" ~ "chemoVehicle_1",
    hash.ID == "TotalSeq-A0302-GGTCGAGAGCATTCA" ~ "chemoVehicle_2",
    hash.ID == "TotalSeq-A0303-CTTGCCGCATGTCAT" ~ "chemoDay2_1",
    hash.ID == "TotalSeq-A0305-CTTTGTCTTTGTGAG" ~ "chemoDay2_2",
    #TRUE ~ hash.ID  # This line means "if none of the above conditions are true, keep the original value"
  )) %>%
  select(Cell.10X.Barcode, referenceID, "hash.ID") %>%
  rename(cellID = Cell.10X.Barcode, barcode = referenceID, sample = "hash.ID")
#write_csv(POOL_1, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_1_map_processed.csv")

POOL_2 = read.csv("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_2_HTO_maxID.csv")
POOL_2_map = read.csv("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_2_map.csv")

POOL_2 = POOL_2 %>%
  rename(Cell.10X.Barcode = CellName) %>%
  inner_join(POOL_2_map, by = "Cell.10X.Barcode") %>%
  filter(!(hash.ID %in% c("Doublet", "Negative"))) %>%
  mutate(hash.ID = case_when(
    hash.ID == "TotalSeq-A0301-ACCCACCAGTAAGAC" ~ "chemoDay5_1",
    hash.ID == "TotalSeq-A0302-GGTCGAGAGCATTCA" ~ "chemoDay5_2",
    hash.ID == "TotalSeq-A0303-CTTGCCGCATGTCAT" ~ "chemoDay7_1",
    hash.ID == "TotalSeq-A0305-CTTTGTCTTTGTGAG" ~ "chemoDay7_2",
    #TRUE ~ hash.ID  # This line means "if none of the above conditions are true, keep the original value"
  )) %>%
  select(Cell.10X.Barcode, referenceID, "hash.ID") %>%
  rename(cellID = Cell.10X.Barcode, barcode = referenceID, sample = "hash.ID")
#write_csv(POOL_2, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_2_map_processed.csv")





################# Formatting map.csv files ###########################################
library(dplyr)
library(readr)
library(stringr)


# Set the working directory to the folder containing your .csv files
setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/retransplant_experiment/")

setwd("C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\SPLINTR\\cell_annotation_tables_cloneInfo\\chemo_experiment\\")

# List all .csv files that have 'map' somewhere in their name
files <- list.files(pattern = "map.*\\.csv$")
files = "C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\SPLINTR\\cell_annotation_tables_cloneInfo\\chemo_experiment\\KRAS_T0_2_map.csv"

# Define a function to process each file
process_file <- function(file_name) {
  # Read the .csv file
  df <- read_csv(file_name)

  # Filter columns
  df <- df %>%
    select(Cell.10X.Barcode, referenceID)

  # Rename columns
  df <- df %>%
    rename(cellID = Cell.10X.Barcode, barcode = referenceID)

  # Extract the sample name from the file name (text before the first "_")
  #sample_name <- str_extract(file_name, "^[^_]+")
  sample_name = "chemoKRASBaseline_2"

  # Create the sample name
  #sample_name <- paste0("inVitro_", sample_name)

  # Add the sample column
  df$sample <- sample_name

  # Write the result back to a csv file
  write_csv(df, str_replace(file_name, ".csv", "_processed.csv"))
}

# Apply the function to each file
#lapply(files, process_file)

process_file("C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\SPLINTR\\cell_annotation_tables_cloneInfo\\chemo_experiment\\KRAS_T0_2_map.csv")
### combining all the files:

library(dplyr)
library(readr)

# Set the main directory where the subfolders are located
main_directory <- "C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\SPLINTR\\cell_annotation_tables_cloneInfo\\"
setwd(main_directory)

# List all subfolders in the main directory
subfolders <- list.dirs(full.names = TRUE, recursive = FALSE)

# Initialize an empty data frame to store the combined data
all_samples <- data.frame()

# Define a function to read and combine csv files from a folder
read_and_combine <- function(folder) {
  # List all _processed.csv files in the subfolder
  files <- list.files(folder, pattern = "_processed\\.csv$", full.names = TRUE)

  # Read all files and combine them into a single data frame
  data <- lapply(files, read_csv) %>% bind_rows()
  return(data)
}

# Apply the function to each subfolder and combine the results
all_samples <- lapply(subfolders, read_and_combine) %>% bind_rows()

all_samples$cellID <- sub("-1$", "", all_samples$cellID)

# Write the combined data frame to a new CSV file in the main folder
write_csv(all_samples, "all_samples_barcode_counts.csv")



###### Figuring out which MLL 10X dataset to use ####################################

library(readr)
library(dplyr)

# File paths
tsv_file_1 <- "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/GSE161676_RAW/GSM4912548_MLL_T0_barcodes.tsv"
tsv_file_2 <- "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/GSE161676_RAW/GSM4912548_MLL_T0_04_barcodes.tsv"
processed_csv_1 <- "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/clonal_competition_experiment/MLL_T0_map_processed.csv"
processed_csv_2 <- "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/clonal_competition_experiment/MLL_T0_map.2_processed.csv"

# Read the data
data_tsv_1 <- read_tsv(tsv_file_1, col_names = FALSE)
data_tsv_2 <- read_tsv(tsv_file_2, col_names = FALSE)
data_csv_1 <- read_csv(processed_csv_1)
data_csv_2 <- read_csv(processed_csv_2)

# Extract values
values_tsv_1 <- data_tsv_1$X1
values_tsv_2 <- data_tsv_2$X1
values_csv_1 <- unique(data_csv_1$cellID)
values_csv_2 <- unique(data_csv_2$cellID)

# Find intersections
overlap_tsv_1_1 <- length(intersect(values_tsv_1, values_csv_1))
overlap_tsv_2_1 <- length(intersect(values_tsv_2, values_csv_1))
overlap_tsv_1_2 <- length(intersect(values_tsv_1, values_csv_2))
overlap_tsv_2_2 <- length(intersect(values_tsv_2, values_csv_2))

# Print the results
cat("Number of overlapping values with first .tsv file, first .csv:", overlap_tsv_1_1, "\n")
cat("Number of overlapping values with second .tsv file, first .csv:", overlap_tsv_2_1, "\n")
cat("Number of overlapping values with first .tsv file, second .csv:", overlap_tsv_1_2, "\n")
cat("Number of overlapping values with second .tsv file, second .csv:", overlap_tsv_2_2, "\n")














######################## Parsing pooled 10X data ####################################

pool1.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_1/", strip.suffix = TRUE)
pool1.mapPath = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_1_map_processed.csv"

pool2.data = Read10X("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/chemo_demultiplexingHashtags/POOL_2/", strip.suffix = TRUE)
pool2.mapPath = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/cell_annotation_tables_cloneInfo/chemo_experiment/POOL_2_map_processed.csv"


pool2 = CreateSeuratObject(pool2.data)

# Load the mapping file
mapping <- read.csv(pool2.mapPath, stringsAsFactors = FALSE)
mapping$cellID <- str_replace(mapping$cellID, "-1$", "")


# Check the head of the mapping dataframe
head(mapping)

# Get the unique sample names
unique_samples <- unique(mapping$sample)

# Initialize an empty list to store the Seurat objects
seurat_list <- list()

# Loop through each sample and create a Seurat object
for (sample in unique_samples) {
  # Get the cell IDs for this sample
  cell_ids <- unique(mapping$cellID[mapping$sample == sample])

  # Subset the Seurat object to include only cells in this sample
  seurat_list[[sample]] <- subset(pool2, cells = cell_ids)
}


############## Saving 10X counts matrices ##############

for (sample_name in names(seurat_list)) {
  seurat_object <- seurat_list[[sample_name]]

  # Ensure the counts are in a sparse matrix format
  counts_matrix <- as(seurat_object@assays$RNA@counts, "sparseMatrix")

  # Define the output path
  output_path <- paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/SPLINTR/10X_HTOdemux/", sample_name, "_10X")

  # Use DropletUtils to write the 10x counts
  DropletUtils::write10xCounts(
    path = output_path,
    x = counts_matrix
  )
}








