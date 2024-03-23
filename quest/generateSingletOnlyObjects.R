### making singlet-only objects for use with doublet functional consequence assessments for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20231220
### Last edited by Madeline E Melzer on 20240213

library(tidyverse)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(glue)
library(umap)
#library(scperturbR)
library(vegan)
library(SingleCellExperiment)
library(ggpubr)
options(future.globals.maxSize = 4000 * 1024^2)

set.seed = 23

srcDirectory = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/"

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


##########################################################################################################################
# making singlet-only objects
##########################################################################################################################


### function for extracting singlets
extractSinglets = function(dataset, sample) {
  rnaDirectory = glue("{srcDirectory}/{dataset}/10X/{sample}/")
  #singletsListDirectory = glue("{srcDirectory}/{dataset}/fatemapID/")
  singlets_file = glue("{rnaDirectory}singlets_all.txt")
  
  expression_matrix <- Read10X(data.dir = rnaDirectory)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  
  singlets <- read.delim(singlets_file, header = FALSE, sep = "\n")[, 1] %>%
    as.vector()
  
  # Loop through the numbers -1 to -6 to check for each suffix
  for (i in 1:6) {
    suffix <- paste0("-", i)  # Create the suffix string ("-1", "-2", ..., "-6")
    
    # Check if the current suffix is part of the 10X barcodes
    if (grepl(suffix, seurat_object |> colnames() |> head(1))) {
      # Read singlets and append the current suffix
      singlets <- read.delim(singlets_file, header = FALSE, sep = "\n")[, 1] |>
        as.vector() |>
        paste0(suffix)
      break  # Exit the loop once a matching suffix is found
    }
  }
  
  # If no matching suffix was found, remove potential suffixes "-1" to "-6"
  if (is.null(singlets)) {
    singlets <- read.delim(singlets_file, header = FALSE, sep = "\n")[, 1] |>
      as.vector()
    for (i in 1:6) {
      singlets <- singlets |> str_remove(paste0("-", i))
    }
  }
  
  #subset Seurat object so it only contains singlets
  singletsSeurat = subset(seurat_object, cells = singlets)
  
  # Count the number of singlets present in the Seurat object
  numSinglets <- length(x = colnames(singletsSeurat))
  
  singletObjectsDirectory = "~/ZhangMelzerEtAl/singletObjects/"
  # Create dataset directory if it does not exist yet
  datasetSingletDirectory <- glue("{singletObjectsDirectory}{dataset}/")
  if (!dir.exists(datasetSingletDirectory)) {
    dir.create(datasetSingletDirectory, recursive = TRUE)
  }
  
  singletObjectPath<- glue("{datasetSingletDirectory}{sample}_singletsOnly.rds")
  # Save the RDS file if it does not exist yet
  if (!file.exists(singletObjectPath)) {
    saveRDS(singletsSeurat, file = singletObjectPath)
  }
  
  return(list(singletsSeurat = singletsSeurat, numSinglets = numSinglets))
}


##########################################################################################################################
# making all the singlet RDS objects and counting singlets per sample
##########################################################################################################################

singletCounts = data.frame(dataset = character(), sample = character(), numSinglets = integer(), stringsAsFactors = FALSE)

datasetDirs = list.dirs("/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data") 
datasetDirs <- grep("^/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/[^/]+/?$", datasetDirs, value = TRUE)
datasetDirs = c("/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/SPLINTR/") 
# Loop through each dataset directory
for (datasetDir in datasetDirs) {
  datasetName = basename(datasetDir)
  rnaDir = file.path(datasetDir, "10X")
  
  if (dir.exists(rnaDir)) {
    sampleDirs = list.dirs(rnaDir, recursive = FALSE, full.names = TRUE)
    sampleDirs = c("inVitro_MLL_1")
    for (sampleDir in sampleDirs) {
      sampleName = basename(sampleDir)
      
      # Extract singlets for current dataset and sample
      singletCount = extractSinglets(datasetName, sampleName)
      
      # Append the dataset, sample, and number of singlets to the results data frame
      singletCounts = rbind(singletCounts, data.frame(dataset = datasetName, sample = sampleName, numSinglets = singletCount$numSinglets))
      
      print(glue("Dataset: {datasetName}, Sample: {sampleName}, Number of Singlets: {singletCount$numSinglets}"))
    }
  }
}

#write.csv(singletCounts, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/singletCountBySample.csv", row.names= FALSE)





##########################################################################################################################
# merge Jain et al. ("Biorxiv") DMSO A and B to make it have over 1000 cells (only 800 cells in smallest sample)
##########################################################################################################################

sample1_dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/Biorxiv/10X/1_DMSO_A"
sample2_dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/Biorxiv/10X/2_DMSO_B"

dataset = "Biorxiv"
sample = "1_DMSO_A"
#setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/fatemapID/")
#sample1 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.01)
sample1 = Read10X(data_dir)
sample1 = CreateSeuratObject(sample1)
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

sample1 = PercentageFeatureSet(sample1, pattern = "^MT-", col.name = "percent.mt")
sample2 = PercentageFeatureSet(sample2, pattern = "^MT-", col.name = "percent.mt")

s1s2_scTransform <- merge(sample1, y = sample2, project = "S1S2")

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "sample")
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

# saving
setwd("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/Biorxiv/10X/DMSO_AB/")
matrix <- s1s2_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1s2_scTransform.integrated@meta.data$label
samples = s1s2_scTransform.integrated@meta.data$sample

Matrix::writeMM(matrix, file = "matrix.mtx")
#system("gzip matrix.mtx")
write.table(features, file = gzfile("features.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile("barcodes.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)





##########################################################################################################################
# troubleshooting LARRY cross-sample cellIDs
##########################################################################################################################


dataset = "LARRY"
sample = "LSK2_d2_3"
larry = read_csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/data/LARRY/LARRY_Barcodes.csv")

larry_filtered = larry %>% filter(sample == "LSK2_d2_3")

length(unique(larry_filtered$cellID))


data_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/10X/", sample, "/")
data_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/LARRY/10X/LSK2_d2_3/corrected_singlet_pairs.csv")

singlets_pair_for_doublet_simulation_file<-list.files(data_dir, pattern="corrected_singlet_pairs.csv")
singlet_pair_df<-read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/LARRY/10X/LSK2_d2_3/corrected_singlet_pairs.csv")


# Extract the sequences or identifiers from columns 1 and 2
column1_sequences <- singlet_pair_df[[1]]  # Replace 1 with the actual column name if it has one
column2_sequences <- singlet_pair_df[[2]]  # Replace 2 with the actual column name if it has one

# Find the common sequences between the two columns
common_sequences <- intersect(column1_sequences, column2_sequences)

# Display the common sequences
print(common_sequences)

# If you want to know the number of common sequences
num_common_sequences <- length(common_sequences)
print(num_common_sequences)



larry_3 = read.csv("/Users/mem3579/Downloads/LARRY_State_3.csv")
larry_2 = read.csv("/Users/mem3579/Downloads/LARRY_State_2.csv")

larry_combined = read.csv("/Users/mem3579/Downloads/LARRY_Barcodes (4).csv")

print(nrow(larry_3))
print(nrow(larry_2))
print(nrow(larry_combined))
print(nrow(larry_3) + nrow(larry_2))

print(length(unique(larry_3$cellID)))
print(length(unique(larry_2$cellID)))
print(length(unique(larry_combined$cellID)))
print(length(unique(larry_3$cellID)) + length(unique(larry_2$cellID)))



print(length(unique(larry_3$barcode)))
print(length(unique(larry_2$barcode)))
print(length(unique(larry_combined$barcode)))
print(length(unique(larry_3$barcode)) + length(unique(larry_2$barcode)))


print(length(unique(larry_3$sample)))
print(length(unique(larry_2$sample)))
print(length(unique(larry_combined$sample)))
print(length(unique(larry_3$sample)) + length(unique(larry_2$sample)))

larry_combined = larry_combined %>% select(-X)

write.csv(larry_combined, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/otherBarcodingDatasets/LARRY/LARRY_barcodes_allStates.csv", row.names = FALSE)

larry_combined = larry_combined[, c("cellID", "barcode", "sample")]




##########################################################################################################################
# correcting singlet numbers
##########################################################################################################################


newSingletDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletStats/scale_factor_0.00003"

singletFiles <- list.files(newSingletDirectory, pattern = "_singlets_stats\\.csv")

data_list <- list()
for (singletFile in singletFiles) {
  df = read_csv(glue("{newSingletDirectory}/{singletFile}"))
  file_name <- basename(singletFile)
  parts <- str_split(file_name, "_", 2)
  dataset <- parts[[1]][1]
  sample_with_suffix <- parts[[1]][2]
  
  # Remove the "_singlets_stats.csv" suffix to get the sample name
  sample <- str_replace(sample_with_suffix, "_singlets_stats\\.csv", "")
  
  df$sample = sample
  data_list[[file_name]] <- df
}


combined_data <- bind_rows(data_list, .id = "file_name")

singletNumbers = select(combined_data, c("dataset", "sample", "total_singlets"))
names(singletNumbers)[names(singletNumbers) == "total_singlets"] <- "numSingletsNEW"
#write.csv(singletNumbers, glue("{newSingletDirectory}/newSingletCountBySample.csv"))

singletNumbers

oldSingletNumbers = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletObjects/singletCountBySample.csv")
names(singletNumbers)[names(singletNumbers) == "numSinglets"] <- "numSingletsOLD"
singletNumbers$dataset <- as.character(singletNumbers$dataset)
singletNumbers$sample <- as.character(singletNumbers$sample)
oldSingletNumbers$dataset <- as.character(oldSingletNumbers$dataset)
oldSingletNumbers$sample <- as.character(oldSingletNumbers$sample)

# Merge the two data frames by 'dataset' and 'sample'
mergedData <- merge(singletNumbers, oldSingletNumbers, by = c("dataset", "sample"), all = TRUE)
mergedData$singletsDiff = mergedData$numSingletsNEW - mergedData$numSingletsOLD
mergedData$percentErrorOnOriginalSinglets = mergedData$singletsDiff / mergedData$numSingletsNEW



#allsampletwobarcode = read.table("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletStats/scale_factor_0.00003/ClonMapper_all_sample_two_barcode_cell_id.txt")

newSingletDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletStats/scale_factor_0.00003"

twoBarcodeCellIDFiles <- list.files(newSingletDirectory, pattern = "_all_sample_two_barcode_cell_id\\.txt")

twoBCList <- list()
for (twoBarcodeCellIDFile in twoBarcodeCellIDFiles) {
  df = read_csv(glue("{newSingletDirectory}/{twoBarcodeCellIDFile}"))
  file_name <- basename(twoBarcodeCellIDFile)
  parts <- str_split(file_name, "_", 2)
  dataset <- parts[[1]][1]
  
  df$dataset = dataset
  twoBCList[[file_name]] <- df
}

twoBCCombined <- bind_rows(twoBCList, .id = "file_name")



##########################################################################################################################
# adding in multi-sample multi-barcode singlets
##########################################################################################################################
newSingletDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletStats/scale_factor_0.00003"

sampleMultiBCSinglets <- list.files(newSingletDirectory, pattern = "_multi_barcode_singlets\\.txt")
datasetMultiBCSinglets <- list.files(newSingletDirectory, pattern = "_all_sample_two_barcode_singlets\\.txt")

### find all unique singlets in the multi-barcode singlet files for a sample



### intersect this with the multi-sample multi-barcode singlet list

### count the number of intersecting cells and add these to the final number of singlets


crossSampleMultiBCSinglets <- numeric()

# Get a list of unique dataset names from your files
#datasets <- unique(sapply(strsplit(list.files(newSingletDirectory, pattern = "__multi_barcode_singlets\\.txt"), "_"), `[`, 1))
datasets = c("Biorxiv", "cellTag","ClonMapper", "FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "LARRY", "non_cancer", "smartseq3_reads", "smartseq3_umis", "SPLINTR", "TREX", "watermelon")

dataset = "ClonMapper"

# Loop through each dataset
for (dataset in datasets) {
  # Get all sampleMultiBCSinglets files for the current dataset
  samplePattern <- paste0("^", dataset, "(_.+)?__multi_barcode_singlets\\.txt$")
  datasetPattern <- paste0("^", dataset, "(_.+)?_all_sample_two_barcode_singlets\\.txt$")
  
  # Get all relevant files using the refined patterns
  sampleFiles <- list.files(newSingletDirectory, pattern = samplePattern)
  print(sampleFiles)
  datasetFile <- list.files(newSingletDirectory, pattern = datasetPattern)
  print(datasetFile)
  
  sampleCellIDs <- unique(unlist(lapply(file.path(newSingletDirectory, sampleFiles), readLines)))
  sampleCellIDsNonUNique = unlist(lapply(file.path(newSingletDirectory, sampleFiles), readLines))
  # Process dataset file
  datasetCellIDs <- readLines(file.path(newSingletDirectory, datasetFile))
  
  # Find overlapping cell IDs
  overlappingIDs <- intersect(sampleCellIDs, datasetCellIDs)
  
  # Calculate and store the final number for the current dataset
  crossSampleMultiBCSinglets[[dataset]] <- length(datasetCellIDs) - length(overlappingIDs)
}

# Print or store the final numbers for each dataset
print(crossSampleMultiBCSinglets)


df_crossSampleMultiBCSinglets <- data.frame(
  dataset = names(crossSampleMultiBCSinglets),
  crossSampleMultiBCCount = unlist(crossSampleMultiBCSinglets),
  stringsAsFactors = FALSE  # to keep dataset names as character strings
)
print(df_crossSampleMultiBCSinglets)
