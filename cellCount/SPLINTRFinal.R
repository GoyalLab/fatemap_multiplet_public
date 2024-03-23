#Importing required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
seed.use  = 10101

#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/SPLINTR/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll <- read.csv(paste0(inputDirectoryPaper, "all_samples_barcode_counts.csv"))

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/SPLINTR/"
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

samples <- c("retransplant","chemoDay2_1", "chemoDay2_2", "chemoDay5_1", "chemoDay5_2", "chemoDay7_1",
             "chemoDay7_2", "chemoKRASBaseline_1", "chemoKRASBaseline_2", "chemoVehicle_1",
             "chemoVehicle_2", "inVitro_FLT3", "inVitro_KRAS", "inVitro_MLL_1", "inVitro_MLL_2")

for (index in 1:length(samples)){
  #Input 10X matrices from paper
  sampleName <- samples[index]
  print(sampleName)
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", sampleName,"/"), strip.suffix = TRUE)
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 3, min.features = 100)
  samplePaper$cellID <- Cells(samplePaper)
  samplePaper[["percent.mito"]] <- PercentageFeatureSet(samplePaper, pattern = "^mt-")
  sampleNameList <- append(sampleNameList, sampleName)
  
  #sample QC filtering
  featuresLower = 100
  featuresUpper = 10000
  countsLower = 1000
  countsUpper = 100000
  percentMito = 20
  
  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, featuresLower, featuresUpper, countsLower, countsUpper, percentMito))
  
  #Getting the list of barcodes for this sample
  barcodeSample <- subset(barcodesAll, sample == sampleName)

  #Creating the standard cell-barcode table for all the datasets
  cellBarcodeDf <- matchBarcodesToCells(samplePaper$cellID, barcodeSample)

  #Filtering out cells which have no associated barcodes
  samplePaper$barcode <- cellBarcodeDf$barcode
  samplePaperBC <- subset(samplePaper, subset = barcode != "not.detected")
  #Data 3
  numBarcodedCellsBeforeUmi <- append(numBarcodedCellsBeforeUmi, as.integer(length(samplePaperBC$cellID)))
  #Data 4
  numBarcodedCellsQCBeforeUmi <- append(numBarcodedCellsQCBeforeUmi, filterQC(samplePaperBC, featuresLower, featuresUpper, countsLower, countsUpper,  percentMito))
  
  umiBarcodedCellDataframe <- umiFilter(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)  
  #Data 5
  numBarcodedCells <- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 6
  numBarcodedCellsQC <- append(numBarcodedCellsQC, filterQC(samplePaper, featuresLower, featuresUpper, countsLower, countsUpper,  percentMito))
  
  #Singlet filtering: Function from the paper
  samplePaper <- findDoubletsByBarcode(samplePaper)
  samplePaper <- subset(samplePaper, doubletBarcode == "singlet")
  #Data 7
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 8
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresLower, featuresUpper, countsLower, countsUpper,  percentMito, save = TRUE, outputPath = paste0(outputRDSPaper, "SPLINTR_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 9
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  sampleSingletCode[["percent.mito"]] <- PercentageFeatureSet(sampleSingletCode, pattern = "^mt-")
  #Data 10
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresLower, featuresUpper, countsLower, countsUpper,  percentMito, save = TRUE, outputPath = paste0(outputRDSSingletCode, "SPLINTR_", sampleName, ".rds")))
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
write.csv(CellNumbers, paste0(outputNumbers, "SPLINTR.csv"))

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresLower, featuresUpper, countsLower, countsUpper,  percentMito, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nFeature_RNA < featuresUpper & nCount_RNA > countsLower & nCount_RNA < countsUpper & percent.mito < percentMito)
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
      df.2 <- data.frame(cell.id = id, barcode = df[1,2], numBarcode = unique.barcodes)
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

# Doublet filtering method from https://atlassian.petermac.org.au/bitbucket/projects/DAW/repos/splintr_paper_code/browse/scRNA_common_scripts/sc-funcs.R
findDoubletsByBarcode <- function(obj, threshold = 2){
  doublets.by.barcode <- c()
  
  # assert barcodes are present in object
  if(is.null(obj$barcode)){
    message("no barcode field identified in object metadata. Stopping")
    break()
  }
  
  # identify cells that have more than 1 barcode
  multbarcodes<-obj$barcode[grep(";",obj$barcode)]
  sortbc <- c()
  for (i in 1:length(multbarcodes)){
    barcodes <- str_split(multbarcodes[i], pattern = ";")
    barcodes <- lapply(barcodes,sort)
    barcodes <- paste(unlist(barcodes),collapse=";")
    names(barcodes) <- names(multbarcodes[i])
    sortbc <- c(sortbc,barcodes)
    
  }
  
  #combination of barcodes that only appear once
  sortbc <- as.data.frame(sortbc)
  doubletbc <- sortbc[which(!(duplicated(sortbc) | duplicated(sortbc, fromLast=TRUE))),,drop=F]
  
  # find barcodes that
  
  #insert metadata into seruat obj
  obj$doubletBarcode <- ifelse(rownames(obj@meta.data) %in% rownames(doubletbc), yes = "doublet", no = "singlet")
  unknown <- which(obj$barcode == "not.detected")
  obj$doubletBarcode[unknown] <- "unknown"
  
  return(obj)
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

######################
#Checking one bug - Error in readRDS for inVitro_MLL_1_singletsOnly.rds
rdsFilePath <- "path/to/your/file.rds"  # Replace with your actual file path
if(file.exists(rdsFilePath)) {
  cat("The file exists. Attempting to read...\n")
  sampleSingletCode <- readRDS(rdsFilePath)
  cat("File read successfully.\n")
} else {
  cat("The file does not exist at the specified path.\n")
}


temp_file_path <- tempfile()
file.copy(rdsFilePath, temp_file_path)
sampleSingletCode <- readRDS(temp_file_path)
