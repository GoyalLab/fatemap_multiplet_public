#Importing required libraries
library(Seurat)
library(stringr)
library(dplyr)
library(tidyverse)
seed.use = 10101


#Input and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/LARRY/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll_1 <- read_csv(paste0(inputDirectoryPaper,"LARRY_State_2.csv"))
barcodesAll_2 <- read_csv(paste0(inputDirectoryPaper,"LARRY_State_3.csv"))
barcodesAll <- rbind(barcodesAll_1, barcodesAll_2)
barcodesAll$cellID <- gsub("-", "", barcodesAll$cellID)

paperFinalCells <- read.table(file = "~/Keerthana/CellTypeCountData/AdditionalData/LARRY_FinalData/stateFate_inVitro_metadata.txt", sep = "\t", header = T)
paperFinalCells$Cell.barcode <- gsub("-", "", paperFinalCells$Cell.barcode)

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/LARRY/"
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
numBarcodedCells <- c()
numBarcodedCellsQC <- c()
numPaperSinglets <- c()
numPaperSingletsQC <- c()

samples <- c("d2_1","d2_2","d2_3","LK_d2", "LK_d4_1", "LK_d4_2", "LK_d6_1_1", "LK_d6_1_2", "LK_d6_2_1",
             "LK_d6_2_2","LSK_d2_1", "LSK_d2_2", "LSK_d2_3", "LSK_d4_1_1", "LSK_d4_1_2", 
             "LSK_d4_1_3", "LSK_d4_2_1", "LSK_d4_2_2", "LSK_d4_2_3", "LSK_d6_1_1", 
             "LSK_d6_1_2", "LSK_d6_1_3", "LSK_d6_2_1", "LSK_d6_2_2", "LSK_d6_2_3")
unmatched <- c("d4_L_1", "d4_L_2", "d4_R_1", "d6_L_1", "d6_L_2","d6_R_1", "d6_R_2", "d6_R_3")


for (index in 1:length(samples)){
  #Input 10X matrices from paper
  sampleName <- samples[index]
  print(sampleName)
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", sampleName,"/"))
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 3, min.features = 200)
  samplePaper$cellID <- Cells(samplePaper)
  samplePaper[["percent.mito"]] <- PercentageFeatureSet(samplePaper, pattern = "^mt-")
  sampleNameList <- append(sampleNameList, sampleName)
  
  #sample QC filtering
  percentMito = 20
  
  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, percentMito))
  
  #Getting the list of barcodes for this sample
  barcodeSample <- subset(barcodesAll, sample == sampleName)
  #Creating the standard cell-barcode table for all the datasets
  cellBarcodeDf <- matchBarcodesToCells(samplePaper$cellID, barcodeSample)
  write.csv(cellBarcodeDf, paste0(outputCellBarcodeDf, "LARRY_", sampleName, "_CellDF.csv" ))
  
  #Filtering out cells which have no associated barcodes
  samplePaper <- subset(samplePaper, subset = cellID %in% barcodeSample$cellID)
  #Data 3
  numBarcodedCells <- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 4
  numBarcodedCellsQC <- append(numBarcodedCellsQC, filterQC(samplePaper, percentMito))
  
  # #Singlet filtering: If cell ID is in their final list, it is a singlet
  paperFinalCellsSample <- subset(paperFinalCells, Library == sampleName)
  cellBarcodeDfSinglets <- as.data.frame(intersect(samplePaper$cellID, paperFinalCellsSample$Cell.barcode))
  samplePaper <- subset(samplePaper, subset = cellID %in% cellBarcodeDfSinglets$`intersect(samplePaper$cellID, paperFinalCellsSample$Cell.barcode)`)

  # # #Data 5
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  # # #Data 6
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, percentMito, save = TRUE, outputPath = paste0(outputRDSPaper, "LARRY_", sampleName, ".rds")))
  
  # #Data 5 and 6 for unmatched is 0 because there is no equivalent sample in main paper (renamed)
  # numPaperSinglets <- append(numPaperSinglets, 0)
  # numPaperSingletsQC <- append(numPaperSingletsQC, 0)
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 7
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  sampleSingletCode[["percent.mito"]] <- PercentageFeatureSet(sampleSingletCode, pattern = "^mt-")
  #Data 8
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, percentMito, save = TRUE, outputPath = paste0(outputRDSSingletCode, "LARRY_", sampleName, ".rds")))
}


CellNumbers <- cbind(tibble(sampleNameList), 
                     tibble(num10xCells), 
                     tibble(num10xCellsQC),
                     tibble(numBarcodedCells), 
                     tibble(numBarcodedCellsQC),
                     tibble(numPaperSinglets),
                     tibble(numPaperSingletsQC),
                     tibble(numTrueSinglets),
                     tibble(numTrueSingletsQC))

write.csv(CellNumbers, paste0(outputNumbers, "LARRY.csv"))

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, percentMito, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = percent.mito < percentMito)
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(nCells)
}

#Method to match barcodes to cells in a Seurat object
matchBarcodesToCells <- function(CellList, BarcodeDf){
  # matches cells in a single cell experiment to detected DNA barcodes.
  # all cells is a list of all cells in the experiment
  # BarcodeDf is a dataframe with cell id and barcode id as columns
  # dataframe returned will have all cells matched to a barcode
  # if there is no barcode matchable to a cell "not.detected" is returned
  # for cells that have multiple detected barcodes each barcode is returned separated by ';'
  cellBarcodeAnnotation <- data.frame()
  for (id in CellList){
    keep <- which(BarcodeDf$cellID == id)
    df <- BarcodeDf[keep,]
    unique.barcodes <- length(unique(df$barcode))
    if (unique.barcodes == 0){
      df.2 <- data.frame(cell.id = id, barcode = "not.detected", numBarcode = unique.barcodes)
      cellBarcodeAnnotation <- rbind(df.2, cellBarcodeAnnotation)
    }
    if (unique.barcodes == 1){
      df.2 <- data.frame(cell.id = id, barcode = df$barcode, numBarcode = unique.barcodes)
      cellBarcodeAnnotation <- rbind(df.2, cellBarcodeAnnotation)
    }
    if (unique.barcodes > 1){
      barcodes <- paste(unique(df$barcode), collapse = ";")
      df.2 <- data.frame(cell.id = id, barcode = barcodes, numBarcode = unique.barcodes)
      cellBarcodeAnnotation <- rbind(df.2, cellBarcodeAnnotation)
    }
  }
  return(cellBarcodeAnnotation)
}

