#Importing required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
seed.use  = 10101

#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/ClonMapper/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll <- read.csv(paste0(inputDirectoryPaper, "all_samples_barcode_counts.csv"))

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/ClonMapper/"
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


samples <- c("TP0", "TP0_5", "FM1", "FM7")
index = 1
for (index in 1:length(samples)){
  #Input 10X matrices from paper
  sampleName <- samples[index]
  print(sampleName)
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", sampleName,"/"))
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 3, min.features = 200)
  samplePaper$cellID <- Cells(samplePaper)
  samplePaper[["percent.mito"]] <- PercentageFeatureSet(samplePaper, pattern = "^MT-")
  sampleNameList <- append(sampleNameList, sampleName)
  
  #sample QC filtering
  featuresLower = 750
  countsLower = 4000
  percentMito = 8
  
  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, featuresLower, countsLower, percentMito))
  
  #Getting the list of barcodes for this sample
  barcodeSample <- subset(barcodesAll, sample == sampleName)
  barcodeSample$cellID <- paste0(barcodeSample$cellID, "-1")

  #Filtering out cells which have no associated barcodes
  samplePaperBC <- subset(samplePaper, subset = cellID %in% barcodeSample$cellID)
  #Data 3
  numBarcodedCellsBeforeUmi <- append(numBarcodedCellsBeforeUmi, as.integer(length(samplePaperBC$cellID)))
  #Data 4
  numBarcodedCellsQCBeforeUmi <- append(numBarcodedCellsQCBeforeUmi, filterQC(samplePaperBC, featuresLower, countsLower, percentMito))
  
  #Filtering out cells based on UMI counts of barcodes
  umiBarcodedCellDataframe <- umiFilter(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)
  
  #Data 5
  numBarcodedCells<- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 6
  numBarcodedCellsQC<- append(numBarcodedCellsQC, filterQC(samplePaper, featuresLower, countsLower, percentMito))
  
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)
  samplePaperBC <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)
  #Singlet filtering: Nuothing is done - so all cells are singlets.
  #Data 7
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 8
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresLower, countsLower, percentMito, save = TRUE, outputPath = paste0(outputRDSPaper, "ClonMapper_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 9
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  sampleSingletCode[["percent.mito"]] <- PercentageFeatureSet(sampleSingletCode, pattern = "^MT-")
  #Data 10
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresLower, countsLower, percentMito, save = TRUE, outputPath = paste0(outputRDSSingletCode, "ClonMapper_", sampleName, ".rds")))
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
write.csv(CellNumbers, paste0(outputNumbers, "ClonMapper.csv"))

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresLower, countsLower, percentMito, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nCount_RNA > countsLower & percent.mito < percentMito)
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(nCells)
}


umiFilter <- function(df){
  colnames(df) <- c("cellID", "BC50StarcodeD8", "sample")
  df <- df %>%
    group_by(cellID, BC50StarcodeD8, sample) %>%
    mutate(nUMI = n()) %>%
    unique() %>%
    ungroup()
  totUMI <- sum(df$nUMI)
  umiCutoff = max(2, 3*totUMI/(4e5))
  dfGood <- df %>%
    filter(nUMI >= umiCutoff)
  return(dfGood)
}
