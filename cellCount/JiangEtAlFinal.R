library(Seurat)
library(tidyverse)
seed.use  = 10101

#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/FateMap/non-cancer/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll <- as.tibble(read.table(paste0(inputDirectoryPaper, "stepFourStarcodeShavedReads50.txt"), header = TRUE))
colnames(barcodesAll) <- c("cellID", "UMI", "barcode", "sample")
barcodesAll <- barcodesAll %>% select(-"UMI")

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/non_cancer/"
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

samples = c("10kbarsiblingA", "10kbarsiblingB")
index =1
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
  featuresUpper = 10000

  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, featuresUpper))
  
  #Getting the list of barcodes for this sample
  sampleNum <- as.character(index)
  barcodeSample <- subset(barcodesAll, sample == sampleNum)
  barcodeSample$cellID <- paste0(barcodeSample$cellID, "-1")
  
  #Filtering out cells which have no associated barcodes
  samplePaperBC <- subset(samplePaper, subset = cellID %in% barcodeSample$cellID)
  #Data 3
  numBarcodedCellsBeforeUmi <- append(numBarcodedCellsBeforeUmi, as.integer(length(samplePaperBC$cellID)))
  #Data 4
  numBarcodedCellsQCBeforeUmi <- append(numBarcodedCellsQCBeforeUmi, filterQC(samplePaperBC, featuresUpper))
  
  #Filtering out cells based on UMI counts of barcodes
  umiBarcodedCellDataframe <- umiFilter(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)  
  
  #Data 5
  numBarcodedCells <- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 6
  numBarcodedCellsQC <- append(numBarcodedCellsQC, filterQC(samplePaper, featuresUpper))
  
  #Singlet filtering: Based on function in the paper
  singletList <- doubletFiltering(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% singletList$cellID)
  
  #Data 7
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 8
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresUpper, save = TRUE, outputPath = paste0(outputRDSPaper, "JiangEtAl_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 9
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  sampleSingletCode[["percent.mito"]] <- PercentageFeatureSet(sampleSingletCode, pattern = "^MT-")
  #Data 10
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresUpper, save = TRUE, outputPath = paste0(outputRDSSingletCode, "JiangEtAl_", sampleName, ".rds")))
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
write.csv(CellNumbers, paste0(outputNumbers, "JiangEtAl.csv"))

#Function to perform QC on a matrix and return the number of cells left and save the object after QC 
filterQC <- function(sample, featuresUpper, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA < featuresUpper)
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(nCells)
}

#Singlet filtering method from the paper
doubletFiltering <- function(barcodes10X){
  # umiCut = 2 # Minimum UMI cutoff (per barcode) for reliable analysis. NOTE that unless this is taking place after unique() subset of barcodes it is not actually an umiCut but rather cutoff for total number of reads (without accounting for PCR duplicates)
  fracumiCut = 0.3
  #to a point, increasing fracumiCut decreases total cells recovered but increases proportion with 1 barcode
  
  linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.
  barcodes10X <-  barcodes10X %>% 
    group_by(cellID, barcode, sample) %>%
    summarise(nUMI = length(sample)) 
  totUMI = sum(barcodes10X$nUMI)
  umiCut = max(2, round(3*totUMI/(4e5)))
  barcodePostUmiCutoff = barcodes10X %>% 
    filter(nUMI >= umiCut) %>% 
    group_by(cellID,sample) %>%
    mutate(fracUMI = nUMI/sum(nUMI)) %>%
    filter(fracUMI >= fracumiCut) %>%
    group_by(cellID, sample) %>%
    mutate(nLineages = length(cellID)) 
  #6952 for umiCut=2+fracumiCut=0.3
  
  ###taking only single cellid-->barcode mappings; removes doublets, ambient RNA (more than 1 barcode per cellID)
  barcodePostUmiLineageCutoff = barcodePostUmiCutoff %>%
    ungroup() %>%
    filter(nLineages <= linCut) %>%
    unique() 
  return(barcodePostUmiLineageCutoff)
}
umiFilter <- function(df){
  colnames(df) <- c("cellID", "BC50StarcodeD8", "sample")
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
