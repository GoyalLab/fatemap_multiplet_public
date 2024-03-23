#Importing required libraries
library(matrixStats)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(stringr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(svglite)

mt_genes <- c(
  "ENSMUSG00000064336", "ENSMUSG00000064337", "ENSMUSG00000064338", "ENSMUSG00000064339",
  "ENSMUSG00000064340", "ENSMUSG00000064341", "ENSMUSG00000064342", "ENSMUSG00000064343",
  "ENSMUSG00000064344", "ENSMUSG00000064345", "ENSMUSG00000064346", "ENSMUSG00000064347",
  "ENSMUSG00000064348", "ENSMUSG00000064349", "ENSMUSG00000064350", "ENSMUSG00000064351",
  "ENSMUSG00000064352", "ENSMUSG00000064353", "ENSMUSG00000064354", "ENSMUSG00000064355",
  "ENSMUSG00000064356", "ENSMUSG00000064357", "ENSMUSG00000064358", "ENSMUSG00000064359",
  "ENSMUSG00000064360", "ENSMUSG00000064361", "ENSMUSG00000065947", "ENSMUSG00000064363",
  "ENSMUSG00000064364", "ENSMUSG00000064365", "ENSMUSG00000064366", "ENSMUSG00000064367",
  "ENSMUSG00000064368", "ENSMUSG00000064369", "ENSMUSG00000064370", "ENSMUSG00000064371",
  "ENSMUSG00000064372"
)

#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/smartseq3_reads/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll <- read.csv(paste0(inputDirectoryPaper, "all_brains_cellID_barcode_read.csv"))

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/smartseq3_reads/"
outputRDSSingletCode <- "~/Keerthana/CellCounts/ProcessedSingletCode"
outputPlotPaper <- "~/Keerthana/CellCounts/QCPlotsPaper/"
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

samples <- c("brain1", "brain2")

mitoPercentageList = c(5,5)
featuresLowerList = c(2500, 2500)
featuresUpperList = c(11000, 10000)
countsUpperList = c(150000, 130000)
index = 1
for (index in 1:length(samples)){
  sampleName = samples[index]
  print(sampleName)
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", sampleName,"/"))
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 0.001*dim(countsPaper)[2], min.features = 200)
  samplePaper$cellID <- Cells(samplePaper)
  sampleNameList <- append(sampleNameList, sampleName)
  
  #sample QC filtering
  mitoPercent = mitoPercentageList[index]
  featuresLower = featuresLowerList[index]
  featuresUpper = featuresUpperList[index]
  countsUpper = countsUpperList[index]
  
  #Data 1
  num10xCells <- append(num10xCells, as.integer(length(samplePaper$cellID)))
  mtGenesObject <- (intersect(Features(samplePaper), mt_genes))
  samplePaper[["percentMito"]] <- PercentageFeatureSet(samplePaper, features = mtGenesObject)
  # QcPlot <- VlnPlot(object = samplePaper, features = c("nFeature_RNA", "nCount_RNA", "percentMito"), ncol = 3, pt.size = 0.01)
  # QcPlot
  # fScatPlot<-FeatureScatter(samplePaper, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_vline(xintercept=130000) + geom_hline(yintercept=11000) + geom_hline(yintercept=1000)
  # fScatPlot
  # ggsave(plot = QcPlot, filename = paste0(outputPlotPaper, "smartSeq3ReadsQC_", sampleName,".png"), width = 6, height = 6)
  # ggsave(plot = fScatPlot, filename = paste0(outputPlotPaper, "smartSeq3ReadsQC2_", sampleName,".png"), width = 6, height = 6)
  #Data 2
  num10xCellsQC <- append(num10xCellsQC, filterQC(samplePaper, featuresLower, featuresUpper, countsUpper, mitoPercent))
  #Getting the list of barcodes for this sample
  barcodeSample <- subset(barcodesAll, sample == sampleName)

  #Filtering out cells which have no associated barcodes
  samplePaperBC <- subset(samplePaper, subset = cellID %in% barcodeSample$cellID)
  
  #Data 3
  numBarcodedCellsBeforeUmi <- append(numBarcodedCellsBeforeUmi, as.integer(length(samplePaperBC$cellID)))
  #Data 4
  numBarcodedCellsQCBeforeUmi <- append(numBarcodedCellsQCBeforeUmi, filterQC(samplePaperBC, featuresLower, featuresUpper, countsUpper, mitoPercent))
  
  #Filtering out cells based on UMI counts of barcodes
  umiBarcodedCellDataframe <- umiFilter(barcodeSample)
  samplePaper <- subset(samplePaper, subset = cellID %in% umiBarcodedCellDataframe$cellID)  
  #Data 5
  numBarcodedCells <- append(numBarcodedCells, as.integer(length(samplePaper$cellID)))
  #Data 6
  numBarcodedCellsQC <- append(numBarcodedCellsQC, filterQC(samplePaper, featuresLower, featuresUpper, countsUpper, mitoPercent))
  
  #No singlet filtering is done.
  #Data 5
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 6
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresLower, featuresUpper, countsUpper, mitoPercent, save = TRUE, outputPath = paste0(outputRDSPaper, "smartSeq3Reads_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 7
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  mtGenesObject <- (intersect(Features(samplePaper), mt_genes))
  sampleSingletCode[["percentMito"]] <- PercentageFeatureSet(sampleSingletCode, features = mtGenesObject)
  #Data 8
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresLower, featuresUpper, countsUpper, mitoPercent, save = TRUE, outputPath = paste0(outputRDSSingletCode, "smartSeq3Reads_", sampleName, ".rds")))
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

write.csv(CellNumbers, paste0(outputNumbers, "smartSeq3Reads.csv"))

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresLower, featuresUpper, countsLower, mitoPercent, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nFeature_RNA < featuresUpper & nCount_RNA < countsUpper & percentMito < mitoPercent)  
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
