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

###function to preprocess the Seurat Object
preprocessSeurat <- function(seuratObject){
  seuratObject <- SCTransform(seuratObject)
  seuratObject <- RunPCA(seuratObject, ndims = 50)
  seuratObject <- RunUMAP(seuratObject, dims = 1:50, seed = set.seed)
  return(seuratObject)
}

### function for extracting singlets from the sheet
extractSinglets <- function(rnaDirectory, seuratObject){
  #Add identity to seurat object to indicate it is a singlet
  singletsFile = glue("{rnaDirectory}singlets_all.txt")
  
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
  return(singlets)
}

#Function to identify the true singlets, barcodeDoublet, non-barcoded cells and umiFiltered cells in the Seurat object for each sample
processSample <- function(srcDirectory, experiment, sampleName){
  rnaDirectory = glue("{srcDirectory}/{sampleName}/")
  # expression_matrix <- Read10X(data.dir = rnaDirectory, strip.suffix = TRUE)
  seuratObject <- readRDS(glue("~/Keerthana/singletCode/plot10xAll/Newplots/rdsFiles/{experiment}_{sampleName}.rds"))
  # Remove any suffix starting with "-"
  new_cell_ids <- gsub("-.*", "", colnames(seuratObject))
  
  # Assign the new cell names without any suffix to the Seurat object
  colnames(seuratObject) <- new_cell_ids
  
  
  seuratObject[["singletStatus"]] <- "notSinglet"
  barcodeDetails <- getBarcodedCells(experiment, sampleName)
  barcodedCells <- data.frame(cellID = unique(barcodeDetails$cellID))
  seuratObject$singletStatus[!Cells(seuratObject) %in% barcodedCells$cellID] <- "notBarcoded"
  
  umiFailedCellList <- umiFailedCells(barcodeDetails)
  validUmiFailedCellList <- unique(umiFailedCellList$cellID[umiFailedCellList$cellID %in% colnames(seuratObject)])
  length(intersect(umiFailedCellList$cellID, Cells(seuratObject)))
  seuratObject$singletStatus[intersect(Cells(seuratObject),validUmiFailedCellList)] <- "umiFailed"
  # Proprocess the seurat object 
  # seuratObject <- preprocessSeurat(seuratObject)
  singlets <- extractSinglets(rnaDirectory, seuratObject)
  validSinglets <- singlets[singlets %in% colnames(seuratObject)]
  seuratObject$singletStatus[Cells(seuratObject) %in% validSinglets] <- "singlet"
  #Data to Plot
  umap_coords <- Embeddings(seuratObject, "umap")
  df <- data.frame(cellID = colnames(seuratObject), umap_1 = umap_coords[, 1],
                   umap_2 = umap_coords[, 2], singletStatus = seuratObject$singletStatus)
  # plotUmapAgain(df, experiment, sampleName)
  singletSummary <- df %>%
    group_by(singletStatus) %>%
    summarise(count = n()) %>%
    spread(key = singletStatus, value = count, fill = 0) # Spread the singletStatus into separate columns
  
  # Now singletSummary should have columns for each singletStatus and their counts
  
  # Create cellCountInfo without subsetting df repeatedly
  cellCountInfo <- data.frame(
    experiment = experiment,
    sampleName = sampleName,
    numCells = length(seuratObject$orig.ident),
    notBarcodedCells = ifelse("notBarcoded" %in% names(singletSummary), singletSummary$notBarcoded, 0),
    numUmiFailedCells = ifelse("umiFailed" %in% names(singletSummary), singletSummary$umiFailed, 0),
    numSinglets = ifelse("singlet" %in% names(singletSummary), singletSummary$singlet, 0),
    numDoublets = ifelse("notSinglet" %in% names(singletSummary), singletSummary$notSinglet, 0),
    stringsAsFactors = FALSE
  )
  # rdsDirectory = "~/Keerthana/plot10xAll/Newplots/rdsFiles/"
  # rdsObjectPath<- glue("{rdsDirectory}{experiment}_{sampleName}.rds")
  # # Save the RDS file if it does not exist yet
  # saveRDS(seuratObject, file = rdsObjectPath)
  return(cellCountInfo)
}

#getting the barcoded for the sample
getBarcodedCells <- function(experiment, sampleName){
  csvPath <- glue("~/Keerthana/singletCode/plot10xAll/Newplots/fatemap_data/{experiment}/fatemapID/stepFourStarcodeShavedReads50.csv")
  barcodesDf <- read_csv(csvPath)
  barcodesDf <- barcodesDf %>% filter(sample == sampleName)
  return(barcodesDf)
}

umiFailedCells <- function(df){
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
  dfBad = data.frame(cellID = setdiff(df$cellID, dfGood$cellID))
  return(dfBad)
}
########################################
#Another plotting function with df and colors as input
plotUmapCount <- function(srcDirectory, experiment, sampleName) {
    rnaDirectory = glue("{srcDirectory}/{sampleName}/")
    seuratObject <- readRDS(glue("~/Keerthana/singletCode/plot10xAll/Newplots/rdsFiles/{experiment}_{sampleName}.rds"))
  
    umap_data <- FetchData(seuratObject, vars = c("umap_1", "umap_2", "nCount_RNA"))
  
    # Arrange the data so that higher nCount_RNA values are last
    # umap_data <- umap_data %>%
    # arrange(nCount_RNA)

    # Now plot using ggplot2, with higher nCount_RNA values plotted on top
    plot <-  ggplot(umap_data, aes(x = umap_1, y = umap_2, color = nCount_RNA)) +
      geom_point(size = 0.5, shape = 16) +
      scale_color_gradient(low = "lightgrey", high = "darkblue") +
      theme_void() +
      theme(legend.position = "right", 
            plot.background = element_blank(), # This sets the plot background to be blank
            panel.background = element_blank(), # Ensure the panel background is blank
            panel.border = element_blank(), # Removes panel border if present
            panel.grid.major = element_blank(), # Removes major grid lines
            panel.grid.minor = element_blank()) # Removes minor grid lines # Adjust alpha for point transparency if desired

  plot  
  # write_csv(umap_data, file = glue("~/Keerthana/plot10xAll/Newplots/countsPlotData/{experiment}_{sampleName}.csv"))
  # ggsave(plot, file = glue("~/Keerthana/plot10xAll/Newplots/countsPlots/{experiment}_{sampleName}.png"), width = 4, height = 4)
}

#########################################
#Another plotting function with df and colors as input
plotUmapAgain <- function(df, experiment, sampleName) {
  umapPngDir <- "~/Keerthana/singletCode/plot10xAll/Newplots/plots/"
  df_category3 <- subset(df, singletStatus == "notBarcoded" | singletStatus == "umiFailed")
  df_rest <- subset(df, singletStatus != "notBarcoded" & singletStatus != "umiFailed")
  rows <- sample(nrow(df_rest))
  df_rest_mixed <- df_rest[rows, ]
  plotDf <- rbind(df_category3, df_rest_mixed)
  coloursGood <- c("singlet" = "#ef8a62", "notSinglet" = "navy" )
  umapGood = ggplot(df_rest, aes(x = umap_1, y = umap_2, color = singletStatus)) +
    geom_point(size = 0.5, shape = 16) +
    scale_color_manual(values = coloursGood) +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umapGood)
  # ggsave(umapGood, file = glue("{umapPngDir}plotsGood/{experiment}_{sampleName}.png"), width = 4, height = 4)
  coloursAll <- c("singlet" = "#ef8a62", "notSinglet" = "navy", "notBarcoded" = "lightgrey", "umiFailed" = "darkgrey")
  umap = ggplot(plotDf, aes(x = umap_1, y = umap_2, color = singletStatus)) +
    geom_point(size = 0.5, shape = 16) +
    scale_color_manual(values = coloursAll) +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umap)
  # 
  # write_csv(df, file = glue("{umapPngDir}plotData/{experiment}_{sampleName}.csv"))
  # ggsave(umap, file = glue("{umapPngDir}plotsAll/{experiment}_{sampleName}.png"), width = 4, height = 4)
}


############################################################################################################
plotUmapDf <- function(df, experiment, sampleName) {
  umapPngDir <- "~/Keerthana/singletCode/plot10xAll/Newplots/plots/plotsRecoloured/"
  df_category3 <- subset(df, singletStatus == "notBarcoded")
  df_category4 <- subset(df,  singletStatus == "umiFailed")
  df_rest <- subset(df, singletStatus != "notBarcoded" & singletStatus != "umiFailed")
  rows <- sample(nrow(df_rest))
  df_rest_mixed <- df_rest[rows, ]
  plotDf <- rbind(df_category3, df_category4, df_rest_mixed)
  coloursAll <- c("singlet" = "#ef8a62", "notSinglet" = "blue", "notBarcoded" = "lightgrey", "umiFailed" = "#264205")
  umap = ggplot(plotDf, aes(x = umap_1, y = umap_2, color = singletStatus)) +
    geom_point(size = 0.5, shape = 16) +
    scale_color_manual(values = coloursAll) +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umap)
  umapTrial = ggplot(df_category3, aes(x = umap_1, y = umap_2)) +  # Removed color mapping from here
    geom_point(color = 'lightgrey', size = 0.5, shape = 16) +  # Set color directly to lightgrey
    geom_point(data = df_category4, aes(x = umap_1, y = umap_2, color = singletStatus), color = "#93A392", size = 0.5, shape = 16) +
    geom_point(data = df_rest_mixed, aes(x = umap_1, y = umap_2, color = singletStatus), size = 0.5, shape = 16) +
    scale_color_manual(values = c("singlet" = "#7B0828", "notSinglet" = "orange")) +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(),  # Sets the plot background to be blank
          panel.background = element_blank(),  # Ensures the panel background is blank
          panel.border = element_blank(),  # Removes panel border if present
          panel.grid.major = element_blank(),  # Removes major grid lines
          panel.grid.minor = element_blank())  # Removes minor grid lines
  
  
  # To display the plot (if running interactively)
  print(umapTrial)
  # write_csv(df, file = glue("{umapPngDir}plotData/{experiment}_{sampleName}.csv"))
  # ggsave(umapTrial, file = glue("{umapPngDir}{experiment}_{sampleName}.png"), width = 4, height = 4)
}

########################################################
singletCounts = data.frame(dataset = character(), sample = character(), numSinglets = integer(), stringsAsFactors = FALSE)
datasetDirs = list.dirs("~/Keerthana/singletCode/plot10xAll/Newplots/fatemap_data") 
# datasetDirs <- grep("^/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/[^/]+?$", datasetDirs, value = TRUE)
# datasetDirs = c("/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/SPLINTR") 
# Loop through each dataset directory

compiledInfo <- data.frame(
  experiment = character(),
  sampleName = character(),
  numCells = integer(),
  barcodedCells = integer(),
  numUmiFailedCells = integer(),
  numSinglets = integer(),
  numDoublets = integer(),
  stringsAsFactors = FALSE # Avoid converting strings to factors
)

for (datasetDir in datasetDirs) {
  datasetName = basename(datasetDir)

  rnaDir = file.path(datasetDir, "10X")
  
  if (dir.exists(rnaDir)) {
    sampleDirs = list.dirs(rnaDir, recursive = FALSE, full.names = TRUE)
    # sampleDirs = c("inVitro_MLL_1")
    for (sampleDir in sampleDirs) {
      sampleName = basename(sampleDir)
      # plotUmapCount(rnaDir,datasetName, sampleName)
      # Extract singlets for current dataset and sample
      temp <- processSample(rnaDir,datasetName, sampleName)
      compiledInfo <- rbind(compiledInfo, temp)
      # Append the dataset, sample, and number of singlets to the results data frame
      # singletCounts = rbind(singletCounts, data.frame(dataset = datasetName, sample = sampleName, numSinglets = singletCount$numSinglets))
      
      print(glue("Dataset: {datasetName}, Sample: {sampleName}"))
    }
  }
}
# write_csv(compiledInfo, "~/Keerthana/plot10xAll/Newplots/finalGoodData.csv")
umapDataDir <- "~/Keerthana/singletCode/plot10xAll/Newplots/plots/plotData/"

umapSvgDir <- "~/Keerthana/singletCode/plot10xAll/umapSvgNew/"




# Plotting singlets, doublets and non-barcoded cells in same plot
csv_files <- list.files(umapDataDir, pattern = "\\.csv$", full.names = TRUE)

for (csv_file in csv_files) {
  # Extract experiment and sampleName from the filename
  filename <- basename(csv_file)
  parts <- strsplit(filename, "_")[[1]]  # Splitting by underscore
  experiment <- parts[1]
  sampleName <- paste(parts[2:length(parts)], collapse = "_")  # Rejoin the remaining parts
  sampleName <- sub("\\.csv$", "", sampleName)  # Removing the .csv extension
  df <- read_csv(csv_file)
  # Call the plotting function
  plotUmapDf(df, experiment, sampleName)
}
# write_csv(compiledInfo, file = "~/Keerthana/plot10xAll/umapFinal/summaryStats.csv")
