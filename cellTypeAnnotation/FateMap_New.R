library(Seurat)
library(reticulate)
library(dplyr)
library(tidyverse)
library(HGNChelper)
scanorama <- import('scanorama')
seed.use = 10101
#Use the rds files from Madeline 
inputDirectory <- "~/Keerthana/CellCounts/SingletCodeData/"
experimentID10X <- c("FM02", "FM03", "FM04", "FM05", "FM06", "FM08")
outputPlotData <- "~/Keerthana/CellTypeCountData/plotData/UmapData/FateMap/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/FateMap/"
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

qc_info <- tibble(
  experimentId = c("FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08"),
  featuresLower = c(200, 200, 200, 200, 200, 200, 200),
  featuresUpper = c(7200, 7200, 7200, 7000, 7200, 5000, 7200),
  percentMito = c(26, 30, 30, 18, 26, 21, 26)
)

for(experimentID in experimentID10X){
  inputDirectoryExperiment <- paste0(inputDirectory, experimentID, "/")
  mergedSeurat <- inputAllRDS(inputDirectoryExperiment)
  qcMetrics <- filter(qc_info, experimentID == experimentId)
  seuratQCed <- filterQC(mergedSeurat, qcMetrics$featuresLower, qcMetrics$featuresUpper, qcMetrics$percentMito )
  integratedSeurat <- integrateSeuratFunc(seuratQCed)
  rm(seuratQCed)
  rm(mergedSeurat)
  gc()
  preprocessedSeurat <- preprocessSeurat(integratedSeurat)
  rm(integratedSeurat)
  gc()
  cellTypeAnnotationDF <- cellTypeAnnotate(experimentID, preprocessedSeurat)
  write.csv(cellTypeAnnotationDF, paste0(outputPlotData,"GoyalEtAl_", experimentID[-2], ".csv"))
}
#Input to all the samples - take all the 10X for each experiment
inputAllRDS <- function(folderPath){
  rdsFileLists <- list.files(path = folderPath, pattern = "\\.rds$", full.names = TRUE)
  seuratObjectList <- list()
  sampleNameList <- list()
  for (filePath in rdsFileLists) {
    data <- (filePath)
    # Assuming 'data' is already a matrix or a Seurat object. If not, adjust accordingly.
    seuratObject <- readRDS(data)
    seuratObject$cellID <- Cells(seuratObject)
    seuratObject[["percent.mito"]] <- PercentageFeatureSet(object = seuratObject, pattern = "^MT-")
    # Add to list
    splitName <- strsplit(tools::file_path_sans_ext(basename(filePath)), "_")[[1]]
    sampleName <- paste(splitName[-length(splitName)], collapse = "_")
    seuratObject$orig.ident <- sampleName
    seuratObject@project.name <- sampleName
    sampleNameList <- append(sampleNameList, sampleName)
    seuratObjectList[[sampleName]] <- seuratObject
  }
  mergedSeuratObject <- merge(x = seuratObjectList[[1]], y = seuratObjectList[-1], add.cell.ids = sampleNameList)
  return(mergedSeuratObject)
}

integrateSeuratFunc <- function(sampleMerged){
  sampleMerged <- seuratQCed
  datasets <- list()
  geneLists <- list()
  groupNames <- unique(sampleMerged$orig.ident)

  for (groupName in groupNames) {
    # Get the expression data for the specified assay
    assayData <- sampleMerged[["RNA"]][groupName]
    print(dim(assayData))
    # Transpose the data matrix (convert to matrix if not already)
    transposedData <- t(as.matrix(assayData))
    print(dim(transposedData))
    groupIndices <- which(Idents(sampleMerged) == groupName)
    # Store results in lists
    datasets[[groupName]] <- transposedData
    geneLists[[groupName]] <- rownames(assayData)
  }
  
  names(datasets) = NULL
  names(geneLists) = NULL
  seed.use <- as.integer(seed.use)
  integratedCorrectedData = scanorama$correct(datasets, geneLists, return_dimred = TRUE, return_dense = TRUE, seed = seed.use, ds_names = groupNames, verbose = TRUE)
  
  correctedScanorama <- t(do.call(rbind, integratedCorrectedData[[2]]))
  colnames(correctedScanorama) <- colnames(sampleMerged)
  rownames(correctedScanorama) <- integratedCorrectedData[[3]]
  correctedScanoramaPca <- t(do.call(rbind, integratedCorrectedData[[1]]))
  colnames(correctedScanoramaPca) <- colnames(sampleMerged)
  rownames(correctedScanoramaPca) <-paste0("PC", 1:100)
  scanoramaAssay <- CreateAssayObject(data = correctedScanorama)
  sampleMerged[["scanorama"]] <- scanoramaAssay
  sampleMerged[["pca_scanorama"]] <- CreateDimReducObject(embeddings = t(correctedScanoramaPca), key = "PC_", assay = "scanorama")
  DefaultAssay(sampleMerged) <- "scanorama"
  return(sampleMerged)
} 

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresLower, featuresUpper, percentMito, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nFeature_RNA < featuresUpper & percent.mito < percentMito)
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(sample)
}

#Preprocessing for Cell Type Annotation
preprocessSeurat <- function(sampleMerged){
  #Pre-processing
  all.genes <- rownames(sampleMerged)
  sampleMerged <- ScaleData(sampleMerged, features = all.genes)
  sampleMerged <- RunPCA(object = sampleMerged, seed.use = seed.use, assay = "scanorama", features= all.genes, reduction.name = "pca_scanorama", verbose = TRUE)
  
  sampleMerged <- FindNeighbors(object=sampleMerged, features = all.genes, dims = 1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
  sampleMerged <- FindClusters(object=sampleMerged,graph.name = "scanorama_snn", seed.use = seed.use, resolution = 0.8)
  sampleMerged <- RunUMAP(object = sampleMerged, reduction = "pca_scanorama", seed.use = seed.use, dims = 1:50, reduction.name = "umap_scanorama")
  return(sampleMerged)
}
#Cell Type annotation

cellTypeAnnotate <- function(experiment, sample){
  if(experiment == "FM04"){
    db = paste0("~/Keerthana/CellTypeCountData/AdditionalData/Watermelon10X_FinalData/Watermelon10X_scTypeMarkers.xlsx")
    tissue = "Breast"
  }
  else{
    db = paste0("~/Keerthana/CellTypeCountData/AdditionalData/FateMap_FinalData/FateMapMarkers.xlsx")
    tissue = "Melanoma"
  }
  gs_list = gene_sets_prepare(db, tissue)
  
  es.max = sctype_score(scRNAseqData = sample[["scanorama"]]$scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = NULL) 
  
  cL_resutls = do.call("rbind", lapply(unique(sample@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sample@meta.data[sample@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sample@meta.data$seurat_clusters==cl)), 10)
  }))
  
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  
  sample@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    sample@meta.data$customclassif[sample@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  SaveSeuratRds(sample, paste0(outputRDS, experimentID, ".rds"))
  umap_coords <- Embeddings(sample, "umap_scanorama")
  df <- data.frame(umap_coords)
  df$cluster <- sample$customclassif
  colnames(df) <- c("umap_1", "umap_2", "cluster")
  df$cellID <- rownames(df)
  return(df)
}
#Save the final finals in the appropriate location
