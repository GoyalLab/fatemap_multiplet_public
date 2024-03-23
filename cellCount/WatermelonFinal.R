#Inuput and output for the paper method
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/Watermelon/"
outputRDSPaper <- "~/Keerthana/CellCounts/ProcessedPaperMethod/"
outputCellBarcodeDf <- "~/Keerthana/CellCounts/cellBarcodeInfo/"
barcodesAll <- read.csv(paste0(inputDirectoryPaper, "Watermelon_barcode_cellID_UMIcounts.csv"))

#Input and output for singletCode
inputDirectorySingletCode = "~/Keerthana/CellCounts/SingletCodeData/watermelon/"
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


samples = c("T47D-lag-1", "T47D-lag-2", "T47D-late-1", "T47D-late-2", "T47D-naive-1", "T47D-naive-2")
sampleNumbers = c("3", "4", "5", "6", "1", "2")

featuresUpperList = c(8000, 8000, 8000, 8000, 7500, 7500)
featuresLowerList = c(2500, 2500, 2500, 2500, 2500, 2500)
mitoPercentageList = c(25, 25, 25, 25, 20, 20)
countsUpperList = c(60000,60000, 60000, 60000, 60000, 60000)
for (index in 1:length(samples)){
  #Input 10X matrices from paper
  sampleName <- samples[index]
  print(sampleName)
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", sampleName,"/"))
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 3, min.features = 200)
  samplePaper$cellID <- Cells(samplePaper)
  samplePaper[["percent.mito"]] <- PercentageFeatureSet(samplePaper, pattern = "^MT")
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
  sampleNum <- sampleNumbers[index]
  barcodeSample$cellID <- paste0(barcodeSample$cellID, "-", sampleNum)
  
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
  
  
  #singlet filtering - minimum UMI per barcode- cell pair = 3 - doublet filtering is done!
  barcodeSample <- barcodeSample  %>% group_by(cellID, barcode, sample) %>% summarise(nUMI = length(cellID))
  singletList <- barcodeSample %>% filter(nUMI >2)
  samplePaper <- subset(samplePaper, subset = cellID %in% singletList$cellID)
  
  #Data 7
  numPaperSinglets <- append(numPaperSinglets, length(samplePaper$cellID))
  #Data 8
  numPaperSingletsQC <- append(numPaperSingletsQC, filterQC(samplePaper, featuresUpper, featuresLower, countsUpper, percentMito, save = TRUE, outputPath = paste0(outputRDSPaper, "Watermelon_", sampleName, ".rds")))
  
  #Singlet Code part of the dataset
  rdsFilePath = paste0(inputDirectorySingletCode, sampleName, "_singletsOnly.rds")
  print(rdsFilePath)
  sampleSingletCode <- readRDS(rdsFilePath)
  sampleSingletCode$cellID <- Cells(sampleSingletCode)
  #Data 9
  numTrueSinglets <- append(numTrueSinglets, length(sampleSingletCode$cellID))
  sampleSingletCode[["percent.mito"]] <- PercentageFeatureSet(sampleSingletCode, pattern = "^MT")
  #Data 10
  numTrueSingletsQC <- append(numTrueSingletsQC, filterQC(sampleSingletCode, featuresUpper, featuresLower, countsUpper, percentMito, save = TRUE, outputPath = paste0(outputRDSSingletCode, "Watermelon_", sampleName, ".rds")))
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
write.csv(CellNumbers, paste0(outputNumbers, "Watermelon.csv"))

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
