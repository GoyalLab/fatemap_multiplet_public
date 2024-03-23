library(Seurat)
library(dplyr)

outputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/FateMap_FinalData/"
experimentIDList = c("FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08")
inputDirectory <- "/projects/p31666/melzer/ZhangMelzerEtAl/singletObjects/"

Singlets <- data.frame(cellID = character(), sampleName = character(), experiment = character(), stringsAsFactors = FALSE)
write.csv(Singlets, paste0(outputDirectory, "singletList.csv"))
for (experiment in experimentIDList){
  cellIdList <- getCellID(experiment)
  Singlets <- rbind.data.frame(Singlets, cellIdList)
}

getCellID <- function(folderName){
  directoryPath <- paste0(inputDirectory, folderName, "/")
  rdsFiles <- list.files(path = directoryPath, pattern = "\\.rds$", full.names = TRUE)
  cellList <- lapply(rdsFiles, readRdsGetCells)
  combinedCells <- bind_rows(cellList)
  combinedCells$Experiment <- folderName
  combinedCells <- as.data.frame(combinedCells)
  return(combinedCells)
}

readRdsGetCells <- function(filePath) {
  splitName <- strsplit(tools::file_path_sans_ext(basename(filePath)), "_")[[1]]
  
  # Remove the last element and recombine the remaining parts
  sampleName <- paste(splitName[-length(splitName)], collapse = "_")
  print(sampleName)
  seuratObject <- readRDS(filePath)
  cellID <- Cells(seuratObject)
  cellIdDf <- data.frame(cellID = cellID)
  cellIdDf$sample <- sampleName
  return(cellIdDf)
}

