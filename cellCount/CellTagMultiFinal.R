#Importing required libraries
library(matrixStats)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(stringr)
library(dplyr)
library(Matrix)
library(readr)

seed.use = 10101
#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/CellTag/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
barcodesAll <- read_csv(paste0(inputDirectoryPaper,"all_samples.csv"))

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/cellTag/"
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

samples = c("B4D21-RNA-r2-1", "d2-RNA-5")

featuresUpperList = c(6000, 6000)
featuresLowerList = c(1000, 2000)
countsUpperList = c(50000,30000)
mitoPercentageList = c(10,10)
index = 1
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
  featuresUpper = featuresUpperList[index]
  featuresLower = featuresLowerList[index]
  countsUpper = countsUpperList[index]
  percentMito = mitoPercentageList[index]
  
  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, featuresUpper, featuresLower, countsUpper, percentMito))
  
  #Getting the list of barcodes for this sample
  barcodeSample <- subset(barcodesAll, sample == sampleName)
  barcodeSample$cellID <- paste0(barcodeSample$cellID, "-1")
  
  #Filtering out cells which have no associated barcodes
  samplePaperBC <- subset(samplePaper, subset = cellID %in% barcodeSample$cellID)
  #Data 3
  numBarcodedCellsBeforeUmi <- append(numBarcodedCellsBeforeUmi, as.integer(length(samplePaperBC$cellID)))
  #Data 4
  numBarcodedCellsQCBeforeUmi <- append(numBarcodedCellsQCBeforeUmi, filterQC(samplePaperBC, featuresUpper, featuresLower, countsUpper, percentMito))
  
  umiBarcodedCellDataframe <- umiFilter(barcodeSample)
  
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)
  #Data 5
  numBarcodedCells <- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 6
  numBarcodedCellsQC <- append(numBarcodedCellsQC, filterQC(samplePaper, featuresUpper, featuresLower, countsUpper, percentMito))
  
  cellBarcodeDf <- umiBarcodedCellDataframe %>%
    group_by(cellID, sample) %>%
    mutate(numBarcode = n())  %>%
    summarise(numBarcode = n_distinct(BC50StarcodeD8), .groups = 'drop') %>%
    ungroup()
  #Singlet filtering: Number of barcodes per cell must be less than 25 and greater than 2
  cellBarcodeDfSinglets <- cellBarcodeDf %>% filter(numBarcode >= 1 & numBarcode <25)
  samplePaper <- subset(samplePaper, subset = cellID %in% cellBarcodeDfSinglets$cellID)
  
  #Data 7
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 8
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresUpper, featuresLower, countsUpper, percentMito, save = TRUE, outputPath = paste0(outputRDSPaper, "CellTagMulti_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 9
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  sampleSingletCode[["percent.mito"]] <- PercentageFeatureSet(sampleSingletCode, pattern = "^mt-")
  #Data 10
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresUpper, featuresLower, countsUpper, percentMito, save = TRUE, outputPath = paste0(outputRDSSingletCode, "CellTagMulti_", sampleName, ".rds")))
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
write.csv(CellNumbers, paste0(outputNumbers, "CellTagMulti.csv"))

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresUpper, featuresLower, countsUpper, percentMito, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nFeature_RNA < featuresUpper & nCount_RNA < countsUpper & percent.mito < percentMito)
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
    filter(nUMI > umiCutoff)
  return(dfGood)
}


