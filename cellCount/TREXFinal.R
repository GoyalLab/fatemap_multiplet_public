#Importing required libraries
library(matrixStats)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(stringr)
library(dplyr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(svglite)


#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/TREX/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll <- read.csv(paste0(inputDirectoryPaper, "all_brains_umi_count_matrix.csv"))

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/TREX/"
outputRDSSingletCode <- "~/Keerthana/CellCounts/ProcessedSingletCode/"

#Common output with cell numbers for the method
outputNumbers <- "~/Keerthana/CellCounts/finalCellNumbers/"

#Numbers to be counted
sampleNameList <- c()
  #For singlet code
  numTrueSingletsQC <- c()
  numTrueSinglets <- c()
  #For paper method
  num10xCells <- c()
  num10xCellsQC <- c()
  numBarcodedCellsBeforeUmi <- c()
  numBarcodedCellsQCBeforeUmi <- c()
  numBarcodedCells <- c()
  numBarcodedCellsQC <- c()
  numPaperSinglets <- c()
  numPaperSingletsQC <- c()

samples = c("brain1", "brain2", "brain3", "brain4")
# brain1_counts <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", "brain1","/"))
# brain1 <- CreateSeuratObject(counts = brain1_counts, min.cells = 3, min.features = 100)
# brain1$cellID <- Cells(brain1)
# brain2_counts <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", "brain2","/"))
# brain2 <- CreateSeuratObject(counts = brain2_counts, min.cells = 3, min.features = 100)
# brain2$cellID <- Cells(brain2)
# brain3_counts <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", "brain3","/"))
# brain3 <- CreateSeuratObject(counts = brain3_counts, min.cells = 3, min.features = 100)
# brain3$cellID <- Cells(brain3)
# brain4_counts <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", "brain4","/"))
# brain4 <- CreateSeuratObject(counts = brain4_counts, min.cells = 3, min.features = 100)
# brain4$cellID <- Cells(brain4)

# # Find cell IDs that match the pattern
# brain1_ctx_list <- grep("10x39", x = colnames(brain1), value = TRUE)
# brain1_ctx = subset(brain1, subset = cellID %in% brain1_ctx_list)
# brain2_ctx_list <- grep("10x42", x = colnames(brain2), value = TRUE)
# brain2_ctx = subset(brain2, subset = cellID %in% brain2_ctx_list)
# brain3_ctx_list <- grep("10x52", x = colnames(brain3), value = TRUE)
# brain3_ctx = subset(brain3, subset = cellID %in% brain3_ctx_list)
# brain4_ctx_list <- grep("10x65", x = colnames(brain4), value = TRUE)
# brain4_ctx = subset(brain4, subset = cellID %in% brain4_ctx_list)
# ctx <- merge(x = brain1_ctx, y = c(brain2_ctx, brain3_ctx, brain4_ctx),add.cell.ids = c("1", "2", "3", "4"), project = "cortex")
index = 1

for (index in 1:length(samples)){
  #Input 10X matrices from paper
  sampleName = samples[index]
  print(sampleName)
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", sampleName,"/"))
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 0.001*dim(countsPaper)[2], min.features = 200)
  samplePaper$cellID <- Cells(samplePaper)
  sampleNameList <- append(sampleNameList, sampleName)
  
  #sample QC filtering
  featuresLower = 500
  featuresUpper = 10000

  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, featuresLower, featuresUpper))
  
  #Getting the list of barcodes for this sample
  barcodeSample <- subset(barcodesAll, sample = as.character(index))
  
  barcodeSample$cellID <- paste0(barcodeSample$cellID, "-1")

  
  #Filtering out cells which have no associated barcodes
  samplePaperBC <- subset(samplePaper, subset = cellID %in% barcodeSample$cellID)
  #Data 3
  numBarcodedCellsBeforeUmi <- append(numBarcodedCellsBeforeUmi, as.integer(length(samplePaperBC$cellID)))
  #Data 4
  numBarcodedCellsQCBeforeUmi <- append(numBarcodedCellsQCBeforeUmi, filterQC(samplePaperBC, featuresLower, featuresUpper))
  
  umiBarcodedCellDataframe <- umiFilter(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)
  #Data 5
  numBarcodedCells <- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 6
  numBarcodedCellsQC <- append(numBarcodedCellsQC, filterQC(samplePaper, featuresLower, featuresUpper))
  
  #Singlet filtering: Number of barcodes per cell must be less than 25 and greater than 2
  cellBarcodeDfSinglets <- doubletFiltering(samplePaper)
  samplePaper <- subset(samplePaper, subset = cellID %in% cellBarcodeDfSinglets$Cell.ID)
  
  #Data 7
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 8
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresLower, featuresUpper, save = TRUE, outputPath = paste0(outputRDSPaper, "TREX_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 9
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  #Data 10
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresLower, featuresUpper, save = TRUE, outputPath = paste0(outputRDSSingletCode, "TREX_", sampleName, ".rds")))
}

CellNumbers <- cbind(tibble(sampleNameList), 
                     tibble(num10xCells), 
                     tibble(num10xCellsQC),
                     tibble(numBarcodedCellsBeforeUmi), 
                     tibble(numBarcodedCellsQCBeforeUmi),
                     tibble(numBarcodedCells), 
                     tibble(numBarcodedCellsQC),
                     tibble(numPaperSinglets),
                     tibble(numPaperSingletsQC),
                     tibble(numTrueSinglets),
                     tibble(numTrueSingletsQC))
write.csv(CellNumbers, paste0(outputNumbers, "TREX.csv"))

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresLower, featuresUpper, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nFeature_RNA < featuresUpper)  
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(nCells)
}

umiFilter <- function(df){
  colnames(df) <- c("cellID","BC50StarcodeD8", "sample")
  df <- df %>%
    group_by(cellID, BC50StarcodeD8, sample) %>%
    mutate(nUMI = n()) %>%
    unique() %>%
    ungroup()
  totUMI <- sum(df$nUMI)
  umiCutoff = max(2, round(3*totUMI/(4e5)))
  dfGood <- df %>%
    filter(nUMI >= umiCutoff)
  return(dfGood)
}

#Doublet Filtering based on exclusive markers
doubletFiltering <- function(samplePaper){
  exclusiveMarkers <- c("Igf2", "Pf4", "Hexb", "Rsph1", "Pdgfra", "Bmp4", "Mog", "Clic6", "Rgs5", "Cldn5", "Reln", "Igfbpl1", "Slc32a1", "Slc17a7", "Aldoc")
  temp <- subset(samplePaper, features = exclusiveMarkers)

  mat = as.data.frame(temp[["RNA"]]$counts)
  result <- t(apply(mat, 1, calculate_threshold_and_binarize))
  # Binarized matrix
  binarized_matrix <- t(sapply(result, function(x) x[[1]]))
  
  # Stats matrix with row names
  stats_matrix <- t(sapply(result, function(x) unlist(x[[2]])))
  rownames(stats_matrix) <- rownames(mat)
  
  #Finding the number of genes expressed
  column_sums <- colSums(binarized_matrix)
  
  # Get cell names where 0, 1 and more than 1 exclusive marker is expressed
  doublets <- data.frame("Cell ID" = colnames(mat) [column_sums > 1])
  singlets_1 <- colnames(mat) [column_sums == 1]
  singlets_2 <- colnames(mat) [column_sums == 0]
  singlets <- data.frame("Cell ID" = c(singlets_1, singlets_2))
  return(singlets)
}

calculate_threshold_and_binarize <- function(row) {
  threshold = 0
  binarized_row <- ifelse(row > threshold, 1, 0)
  quartile_stats <- c(min = min(row), Q1 = quantile(row, probs = 0.25), median = median(row), mean = mean(row), Q3 = quantile(row, probs = 0.75), max = max(row), threshold = threshold)
  return(list(binarized_row, quartile_stats))
}

