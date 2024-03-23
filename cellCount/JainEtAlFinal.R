library(Seurat)
library(tidyverse)
seed.use  = 10101

#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/FateMap/NJ01/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll <- as_tibble(read.table("~/Keerthana/CellCounts/UnprocessedData/FateMap/NJ01/stepFourStarcodeShavedReads50.txt", header = TRUE))
colnames(barcodesAll) <- c("cellID", "barcode", "sample")

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/Biorxiv/"
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


fileNamePaper = c("filtered_feature_bc_matrix", "filtered_feature_bc_matrix2", 
            "filtered_feature_bc_matrix3", "filtered_feature_bc_matrix4",
            "filtered_feature_bc_matrix5", "filtered_feature_bc_matrix6")
sampleNumList = c("1", "2", "5", "6", "3", "4")
samples <- c("DMSO_A", "DMSO_B", "LSD1i_A", "LSD1i_B", "DOT1Li_A", "DOT1Li_B")

index = 5
for (index in 1:length(samples)){
  #Input 10X matrices from paper
  sampleName <- samples[index]
  print(sampleName)
  matrixName <- fileNamePaper[index]
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", matrixName,"/"))
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 3, min.features = 200)
  samplePaper$cellID <- Cells(samplePaper)
  samplePaper[["percent.mito"]] <- PercentageFeatureSet(samplePaper, pattern = "^MT-")
  sampleNameList <- append(sampleNameList, sampleName)
  
  #sample QC filtering
  featuresLower = 200
  percentMito = 30
  
  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, featuresLower, percentMito))
  
  #Getting the list of barcodes for this sample
  sampleNum <- sampleNumList[index]
  barcodeSample <- subset(barcodesAll, sample == sampleNum)
  barcodeSample$cellID <- paste0(barcodeSample$cellID, "-1")

  #Filtering out cells which have no associated barcodes
  samplePaperBC <- subset(samplePaper, subset = cellID %in% barcodeSample$cellID)
  #Data 3
  numBarcodedCellsBeforeUmi <- append(numBarcodedCellsBeforeUmi, as.integer(length(samplePaperBC$cellID)))
  #Data 4
  numBarcodedCellsQCBeforeUmi <- append(numBarcodedCellsQCBeforeUmi, filterQC(samplePaperBC, featuresLower, percentMito))
  
  #Filtering out cells based on UMI counts of barcodes
  umiBarcodedCellDataframe <- umiFilter(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)  
  
  #Data 5
  numBarcodedCells <- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 6
  numBarcodedCellsQC <- append(numBarcodedCellsQC, filterQC(samplePaper, featuresLower, percentMito))
  
  #Singlet filtering: Based on function in the paper
  singletList <- doubletFiltering(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% singletList$cellID)
  
  #Data 7
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 8
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresLower, percentMito, save = TRUE, outputPath = paste0(outputRDSPaper, "JainEtAl_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, as.character(index), "_", sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 9
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  sampleSingletCode[["percent.mito"]] <- PercentageFeatureSet(sampleSingletCode, pattern = "^MT-")
  #Data 10
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresLower, percentMito, save = TRUE, outputPath = paste0(outputRDSSingletCode, "JainEtAl_", sampleName, ".rds")))
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
write.csv(CellNumbers, paste0(outputNumbers, "JainEtAl.csv"))

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresLower, percentMito, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & percent.mito < percentMito)
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(nCells)
}

#Singlet filtering method from the paper
doubletFiltering <- function(barcodes10X){
  #Checking for pattern for 10X barcodes - ideally all of the barcodes in this file must be 10X - so a check
  barcodes10X %>% .$cellID %>% unique() %>% length()
  barcodes10XFilter <- barcodes10X[grepl("([AT][CG][ATCG]){8,}", barcodes10X$barcode), ]
  barcodes10XFilter %>% .$cellID %>% unique() %>% length()
  
  #Filtering barcodes based on UMI counts and percentage of UMI counts from a barcode vs total (< 0.4 means it is removed)
  fracUMICut = 0.4
  barcodes10XFilter <- barcodes10XFilter %>% group_by(cellID, barcode, sample) %>% summarise(nUMI = length(cellID))
  totUMI <- sum(barcodes10XFilter$nUMI)
  umiCut = max(2, round(3*totUMI/(4e5)))
  upToLineageCounting <- barcodes10XFilter %>% filter(nUMI >= umiCut) %>%
    group_by(cellID, sample) %>% mutate(fracUMI = nUMI/sum(nUMI)) %>% filter(fracUMI > fracUMICut) %>%
    group_by(cellID) %>% mutate(nLineages = length(cellID))
  
  upToLineageCounting %>% .$cellID %>% unique() %>% length()
  
  #Removing cells containing more than one barcode even after the above filtering
  linCut = 1
  linCountToOverlaps <- upToLineageCounting %>% ungroup() %>% filter(nLineages == linCut) %>% unique()
  linCountToOverlaps %>% .$cellID %>% length()
  return(linCountToOverlaps)
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
