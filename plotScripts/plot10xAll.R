remotes::install_version("matrixStats", version="1.1.0") 
library(tidyverse)
library(matrixStats)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(umap)
library(glue)
library(ggplot2)
library(svglite)
library(ggrepel)
library(SingleCellExperiment)
options(future.globals.maxSize = 4000 * 1024^2)
set.seed = 101010

# Identifying singlet only cells in the 10X object, pre-processing it, plotting the UMAP and saving the UMAP files to plot

#Functions
#Function to write umap coordinates along with singlet status
saveUmap <- function(df, experiment, sampleName){
  umapDataDir <- "~/Keerthana/plot10xAll/umapData/"
  write_csv(df, glue("{umapDataDir}{experiment}_{sampleName}.csv"))
}

# Function to plot a UMAP from the seurat object saving the plot to a specified directory and also associated data
plotUmap <- function(seuratObject, experiment, sampleName) {
  umapPngDir <- "~/Keerthana/plot10xAll/umapPNG/"
  umapSvgDir <- "~/Keerthana/plot10xAll/umapSVG/"
  umap_coords <- Embeddings(seuratObject, "umap")
  df <- data.frame(cellID = colnames(seuratObject), umap_1 = umap_coords[, 1],
                      umap_2 = umap_coords[, 2], singletStatus = seuratObject$singletStatus)
  saveUmap(df, experiment, sampleName)
  umapDataDir <- "~/Keerthana/plot10xAll/umapData/"
  df <- read_csv(glue("{umapDataDir}{experiment}_{sampleName}.csv"))
  colours <- c("lightgrey", "hotpink3")
  umap = ggplot(df, aes(x = umap_1, y = umap_2, color = singletStatus)) +
    geom_point(size = 1, shape = 16) +
    scale_color_manual(values = colours) +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umap)
  
  # Save the plot
  svglite(filename = glue("{umapSvgDir}{experiment}_{sampleName}.svg"), width = 4, height = 6)
  plot(umap)
  dev.off()
  ggsave(umap, file = glue("{umapPngDir}{experiment}_{sampleName}.png"))
}

###function to preprocess the Seurat Object
preprocessSeurat <- function(seuratObject){
  seuratObject <- SCTransform(seuratObject)
  seuratObject <- RunPCA(seuratObject, ndims = 50)
  seuratObject <- RunUMAP(seuratObject, dims = 1:50, seed = set.seed)
}

### function for extracting singlets and adding that to the metadata
extractSinglets = function(srcDirectory, dataset, sampleName) {
  rnaDirectory = glue("{srcDirectory}/{sampleName}/")
  singletsFile = glue("{rnaDirectory}singlets_all.txt")
  expression_matrix <- Read10X(data.dir = rnaDirectory, strip.suffix = TRUE)
  seuratObject <- CreateSeuratObject(counts = expression_matrix, min.features = 100)
  
  singlets <- read.delim(singletsFile, header = FALSE, sep = "\n")[, 1] %>%
    as.vector()
  
  # Loop through the numbers -1 to -6 to check for each suffix
  for (i in 1:6) {
    suffix <- paste0("-", i)  # Create the suffix string ("-1", "-2", ..., "-6")
    
    # Check if the current suffix is part of the 10X barcodes
    if (grepl(suffix, seuratObject |> colnames() |> head(1))) {
      # Read singlets and append the current suffix
      singlets <- read.delim(singletsFile, header = FALSE, sep = "\n")[, 1] |>
        as.vector() |>
        paste0(suffix)
      break  # Exit the loop once a matching suffix is found
    }
  }
  
  # If no matching suffix was found, remove potential suffixes "-1" to "-6"
  if (is.null(singlets)) {
    singlets <- read.delim(singletsFile, header = FALSE, sep = "\n")[, 1] |>
      as.vector()
    for (i in 1:6) {
      singlets <- singlets |> str_remove(paste0("-", i))
    }
  }
  
  #Add identity to seurat object to indicate it is a singlet
  valid_singlets <- singlets[singlets %in% colnames(seuratObject)]
  
  seuratObject[["singletStatus"]] <- "not_singlet"
  seuratObject$singletStatus[Cells(seuratObject) %in% valid_singlets] <- "singlet"
  
  # Proprocess the seurat object 
  seuratObject <- preprocessSeurat(seuratObject)
  plotUmap(seuratObject, dataset, sampleName)
  rdsDirectory = "~/Keerthana/plot10xAll/finalRdsFiles/"
  rdsObjectPath<- glue("{rdsDirectory}{dataset}_{sampleName}.rds")
  # Save the RDS file if it does not exist yet
  if (!file.exists(rdsObjectPath)) {
    saveRDS(seuratObject, file = rdsObjectPath)
  }
  # return(list(singletsSeurat = singletsSeurat, numSinglets = numSinglets))
}

############################################
#Function to mark doublets
markDoublet <- function(csvFile){
  filename <- basename(csvFile)
  parts <- strsplit(filename, "_")[[1]]  # Splitting by underscore
  experiment <- parts[1]
  sampleName <- paste(parts[2:length(parts)], collapse = "_")  # Rejoin the remaining parts
  if(experiment == "non"){
      experiment = "non_cancer"
      sampleName = paste(parts[3:length(parts)], collapse = "_")
  }
  if(experiment == "smartseq3"){
    experiment = paste(parts[1:2], collapse = "_")
    sampleName = paste(parts[3:length(parts)], collapse = "_")
  }
  sampleName <- sub("\\.csv$", "", sampleName) 
  df <- read_csv(glue("{umapDataDir}{experiment}_{sampleName}.csv"))
  txt_path <- glue("/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/{experiment}/fatemapID/stepFourStarcodeShavedReads50.txt")
  csv_path <- glue("/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/{experiment}/fatemapID/stepFourStarcodeShavedReads50.csv")
  print(sampleName)
  if(experiment == "TREX"){
    barcodes_df <- read_csv(csv_path)
    barcodes_df$sample = paste0("brain",barcodes_df$sample)
    allCellID <- barcodes_df %>%
      distinct(cellID, sample, .keep_all = TRUE)
    allCellID <- allCellID %>% filter(sample == sampleName)
  }else if (file.exists(txt_path)) {
    barcodes_df <- read_table(txt_path)
    allCellID <- barcodes_df %>%
      distinct(cellID, .keep_all = TRUE)
  } else if (file.exists(csv_path)) {
    barcodes_df <- read_csv(csv_path)
    allCellID <- barcodes_df %>%
    distinct(cellID, sample, .keep_all = TRUE)
    allCellID <- allCellID %>% filter(sample == sampleName)
  } else {
    stop("Neither file exists.")
  }
  if (experiment == "watermelon"){
    df$cellID <- sapply(strsplit(df$cellID, "-"), `[`, 1)
  }
  
  # Assuming 'singletStatus' is a factor or character. If it's not, adjust accordingly.
  
  # Find cell IDs in both dataframes but not marked as singlets
  # Identify cell IDs marked as "singlet" in df
  singlet_ids <- df$cellID[df$singletStatus == "singlet"]
  
  # Find cell IDs that are both in allCellID and df, but not in singlet_ids
  doublet_ids <- allCellID$cellID[allCellID$cellID %in% df$cellID & !(allCellID$cellID %in% singlet_ids)]
  

  # Update the singletStatus to 'doublet' for these IDs
  # Saving both singlet and doulet status
  df$singletStatus[df$cellID %in% doublet_ids] <- "doublet"
  coloursAll <- c("not_singlet" = "lightgrey", "doublet" = "navyblue", "singlet" = "#ef8a62")
  # if(length(singlet_ids) < length(doublet_ids)){
  #   bothDf$singletStatus <- factor(bothDf$singletStatus, levels = c("not_singlet", "doublet", "singlet"))
  # }else{
  #   bothDf$singletStatus <- factor(bothDf$singletStatus, levels = c("not_singlet", "singlet", "doublet"))
  # }
  # bothDf <- bothDf[order(bothDf$singletStatus),]
  plotData <- plotUmapAgain(df, experiment, sampleName, coloursAll, "allThree")
  numbers <- data.frame(
    experiment = experiment,
    sampleName = sampleName,
    barcodedCells = length(intersect(allCellID$cellID, df$cellID)),
    numSinglets = length(singlet_ids),
    numDoublets = length(doublet_ids),
    stringsAsFactors = FALSE
  )
  return(numbers)
  }
#########################################
#Another plotting function with df and colors as input
plotUmapAgain <- function(df, experiment, sampleName, colours, prefix) {
  df_category3 <- df[df$singletStatus == "not_singlet", ]
  df_rest <- df[df$singletStatus != "not_singlet", ]
  rows <- sample(nrow(df_rest))
  df_rest_mixed <- df_rest[rows, ]
  plotDf <- rbind(df_category3, df_rest_mixed)
  umap = ggplot(plotDf, aes(x = umap_1, y = umap_2, color = singletStatus)) +
    geom_point(size = 0.5, shape = 16) +
    scale_color_manual(values = colours) +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umap)
  
  # # # Save the plot
  # svglite(filename = glue("{umapSvgDir}{experiment}_{sampleName}.svg"), width = 4, height = 6)
  # plot(umap)
  # dev.off()
  
  write_csv(df, path = glue("~/Keerthana/plot10xAll/umapFinal/dataUmap/{experiment}_{sampleName}.csv"))
  ggsave(umap, file = glue("{umapPngDir}/{experiment}_{sampleName}.png"), width = 4, height = 4)
  return(plotDf)
}
  

##########################################
#function to just the umap data
plotUmapDf <- function(experiment, sampleName) {
  df <- read_csv(glue("{umapDataDir}{experiment}_{sampleName}.csv"))
  colours <- c("lightgrey", "hotpink3")
  umap = ggplot(df, aes(x = umap_1, y = umap_2, color = singletStatus)) +
    geom_point(size = 0.5, shape = 16) +
    scale_color_manual(values = colours) +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umap)
  
  # # Save the plot
  svglite(filename = glue("{umapSvgDir}{experiment}_{sampleName}.svg"), width = 4, height = 6)
  plot(umap)
  dev.off()
  ggsave(umap, file = glue("{umapPngDir}{experiment}_{sampleName}.png"), width = 4, height = 6)
}

############################################################################################################

singletCounts = data.frame(dataset = character(), sample = character(), numSinglets = integer(), stringsAsFactors = FALSE)
datasetDirs = list.dirs("/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data") 
datasetDirs <- grep("^/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/[^/]+?$", datasetDirs, value = TRUE)
datasetDirs = c("/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/SPLINTR") 
# Loop through each dataset directory
for (datasetDir in datasetDirs) {
  datasetName = basename(datasetDir)
  rnaDir = file.path(datasetDir, "10X")
  
  if (dir.exists(rnaDir)) {
    sampleDirs = list.dirs(rnaDir, recursive = FALSE, full.names = TRUE)
    # sampleDirs = c("inVitro_MLL_1")
    for (sampleDir in sampleDirs) {
      sampleName = basename(sampleDir)
      
      # Extract singlets for current dataset and sample
      extractSinglets(rnaDir,datasetName, sampleName)
      
      # Append the dataset, sample, and number of singlets to the results data frame
      # singletCounts = rbind(singletCounts, data.frame(dataset = datasetName, sample = sampleName, numSinglets = singletCount$numSinglets))
      
      print(glue("Dataset: {datasetName}, Sample: {sampleName}"))
    }
  }
}
umapDataDir <- "~/Keerthana/plot10xAll/umapData/"
umapPngDir <- "~/Keerthana/plot10xAll/umapFinal/plotPngUmap/"
umapSvgDir <- "~/Keerthana/plot10xAll/umapSvgNew/"

compiledInfo <- data.frame(
  experiment = character(),
  sampleName = character(),
  barcodedCells = integer(),
  numSinglets = integer(),
  numDoublets = integer(),
  stringsAsFactors = FALSE # Avoid converting strings to factors
)

# Plotting singlets, doublets and non-barcoded cells in same plot
csv_files <- list.files(umapDataDir, pattern = "\\.csv$", full.names = TRUE)
flag = 1
for (csv_file in csv_files){
  if(flag != 0){
  temp <- markDoublet(csv_file)
  compiledInfo <- rbind(compiledInfo, temp)
  }
}
for (csv_file in csv_files) {
  # Extract experiment and sampleName from the filename
  filename <- basename(csv_file)
  parts <- strsplit(filename, "_")[[1]]  # Splitting by underscore
  experiment <- parts[1]
  sampleName <- paste(parts[2:length(parts)], collapse = "_")  # Rejoin the remaining parts
  sampleName <- sub("\\.csv$", "", sampleName)  # Removing the .csv extension
  
  # Call the plotting function
  plotUmapDf(experiment, sampleName)
}
write_csv(compiledInfo, file = "~/Keerthana/plot10xAll/umapFinal/summaryStats.csv")
