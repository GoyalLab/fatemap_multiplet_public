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
seed.use = 10101

# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

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
inputDirectoryPaper = "~/Keerthana/CellCounts/UnprocessedData/smartseq3_umis/"
outputDirectory <- "~/Keerthana/CellTypeCountData/smartSeq3umi/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/smartSeq3umi/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/UmapData/smartSeq3umi/"

samples <- c("brain1", "brain2")

mitoPercentageList = c(5,5)
featuresLowerList = c(2500, 2500)
featuresUpperList = c(11000, 11000)
countsUpperList = c(150000, 130000)

sampleList <- list()
sampleNameList <- c()
for(index in 1:length(samples)){
  sampleName = samples[index]
  print(sampleName)
  countsPaper <- Read10X(data.dir = paste0(inputDirectoryPaper, "10X/", sampleName,"/"))
  samplePaper <- CreateSeuratObject(counts = countsPaper, min.cells = 0.001*dim(countsPaper)[2], min.features = 200)
  samplePaper$cellID <- Cells(samplePaper)
  mitoPercent = mitoPercentageList[index]
  featuresLower = featuresLowerList[index]
  featuresUpper = featuresUpperList[index]
  countsUpper = countsUpperList[index]
  mitoGenes <- intersect(Features(samplePaper), mt_genes)
  samplePaper[["percentMito"]] <- PercentageFeatureSet(samplePaper, features = mitoGenes)
  samplePaper <- filterQC(samplePaper, featuresLower, featuresUpper, countsUpper, mitoPercent)
  samplePaper@project.name <- sampleName
  sampleNameList <- append(sampleNameList, sampleName)
  sampleList[[sampleName]] <- samplePaper
}
sampleMerged <- merge(x = sampleList[[1]], y = sampleList[-1], add.cell.ids = sampleNameList)
sampleMerged <- NormalizeData(sampleMerged)
sampleMerged <- FindVariableFeatures(sampleMerged, nfeatures = 5000)
sampleMerged <- ScaleData(sampleMerged)
sampleMerged <- RunPCA(object = sampleMerged, npcs = 30, seed.use = seed.use,verbose = TRUE)


sampleMerged <- FindNeighbors(object=sampleMerged, dims = 1:30, force.recalc = TRUE)
sampleMerged <- FindClusters(object=sampleMerged, seed.use = seed.use, resolution = 0.8)
sampleMerged <- RunUMAP(object = sampleMerged,seed.use = seed.use, dims = 1:30)
df <- cellTypeAnnotate(sampleMerged)
write_csv(df, paste0(outputDirectoryData, "smartSeq3umi.csv") )

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresLower, featuresUpper, countsLower, mitoPercent, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nFeature_RNA < featuresUpper & nCount_RNA < countsUpper & percentMito < mitoPercent)  
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(sample)
}

cellTypeAnnotate <- function(sample){
  db = paste0("~/Keerthana/CellTypeCountData/AdditionalData/ModifiedSmartSeq3Markers.xlsx")
  tissue = "Neural"
  gs_list = gene_sets_prepare_local(db, tissue)
  
  es.max = sctype_score(scRNAseqData = sample[["RNA"]]$scale.data, scaled = TRUE, normalized = TRUE,
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
  DimPlot(sample,group.by = "customclassif")
  SaveSeuratRds(sample, paste0(outputRDS,  "smartSeq3umi.rds"))
  umap_coords <- Embeddings(sample, "umap")
  df <- data.frame(umap_coords)
  df$cluster <- sample$customclassif
  colnames(df) <- c("umap_1", "umap_2", "cluster")
  df$cellID <- rownames(df)
  return(df)
}

gene_sets_prepare_local <- function(path_to_db_file, cell_type) {
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,]
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # Clean gene symbols (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i) {
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i], ",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(unique(markers_all))
    if(length(markers_all) > 0) {
      paste0(markers_all, collapse = ",")
    } else {
      ""
    }
  })
  
  # Clean gene symbols (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i) {
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i], ",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(unique(markers_all))
    if(length(markers_all) > 0) {
      paste0(markers_all, collapse = ",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///", ",", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore1 = gsub(" ", "", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///", ",", cell_markers$geneSymbolmore2)
  cell_markers$geneSymbolmore2 = gsub(" ", "", cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]), ","))))
  names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]), ","))))
  names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}
