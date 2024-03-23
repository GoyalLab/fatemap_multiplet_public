lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(svglite)
#remotes::install_version("matrixStats", version="1.1.0")
library(Matrix)
library(ggplot2)
library(shadowtext)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db = paste0("~/Keerthana/CellTypeCountData/AdditionalData/LARRY_FinalData/scTypeMarker-Larry.xlsx")
tissue = "Blood"
gs_list = gene_sets_prepare(db, tissue)

featuresUpperList = c(6000, 6000)
featuresLowerList = c(1000, 2000)
countsUpperList = c(50000,30000)
mitoPercentageList = c(10,10)

#HSC - d2, d5
inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/CellTag_FinalData/matrix_files/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/CellTag/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/CellTag/"
outputDirectory <- "~/Keerthana/CellTypeCountData/CellTagMulti/HSC/"
outputRDS  <- "~/Keerthana/CellTypeCountData/finalRDS/CellTag/"

countsDir_1 <- paste0(inputDirectory, "d2-RNA-5/")
countsDir_2 <- paste0(inputDirectory, "d5-RNA-1/")
countsDir_3 <- paste0(inputDirectory, "d5-RNA-2/")

mat_1 <- readMM(paste0(countsDir_1, "matrix.mtx.gz"))
cells_1 <- read.csv(paste0(countsDir_1, "barcodes.tsv"), sep = "\t", header = F)
features_1 <- read.csv(paste0(countsDir_1, "genes.tsv"), sep = "\t", header = F)
rownames(mat_1) <- features_1$V1
colnames(mat_1) <- cells_1$V1
sample_1 <- CreateSeuratObject(counts = mat_1, min.cells = 3, min.features = 500)

mat_2 <- readMM(paste0(countsDir_2, "matrix.mtx.gz"))
cells_2 <- read.csv(paste0(countsDir_2, "barcodes.tsv"), sep = "\t", header = F)
features_2 <- read.csv(paste0(countsDir_2, "genes.tsv"), sep = "\t", header = F)
rownames(mat_2) <- features_2$V1
colnames(mat_2) <- cells_2$V1
sample_2 <- CreateSeuratObject(counts = mat_2, min.cells = 3, min.features = 500)

mat_3 <- readMM(paste0(countsDir_3, "matrix.mtx.gz"))
cells_3 <- read.csv(paste0(countsDir_3, "barcodes.tsv"), sep = "\t", header = F)
features_3 <- read.csv(paste0(countsDir_3, "genes.tsv"), sep = "\t", header = F)
rownames(mat_3) <- features_3$V1
colnames(mat_3) <- cells_3$V1
sample_3 <- CreateSeuratObject(counts = mat_3, min.cells = 3, min.features = 500)

#QC filtration
sample_1$cellID = Cells(sample_1)
sample_2$cellID = Cells(sample_2)
sample_3$cellID = Cells(sample_3)
sample_1[["percent.mito"]] <- PercentageFeatureSet(sample_1, pattern = "^mt-")
sample_2[["percent.mito"]] <- PercentageFeatureSet(sample_2, pattern = "^mt-")
sample_3[["percent.mito"]] <- PercentageFeatureSet(sample_3, pattern = "^mt-")


sample_1 <- filterQC(sample_1, 6000, 2000, 30000, 10)
sample_2 <- filterQC(sample_2, 6000, 2000, 30000, 10)
sample_3 <- filterQC(sample_3, 6000, 2000, 30000, 10)

#doublet filtration
singletList <- read.csv("~/Keerthana/CellTypeCountData/AdditionalData/CellTag_FinalData/ProcessedData/lsk_rna.csv")
singletList <- singletList %>% mutate(cellID = paste0(cellID, "-1"))
singletCells_1 = data.frame("cellID" = intersect(Cells(sample_1), singletList$cellID))
singletCells_2 = data.frame("cellID" = intersect(Cells(sample_2), singletList$cellID))
singletCells_3 = data.frame("cellID" = intersect(Cells(sample_3), singletList$cellID))



multiplets_1 <- sample_1$cellID[!(sample_1$cellID %in% singletCells_1$cellID)]
multiplets_2 <- sample_2$cellID[!(sample_2$cellID %in% singletCells_2$cellID)]
multiplets_3 <- sample_3$cellID[!(sample_3$cellID %in% singletCells_3$cellID)]

length(unique(sample_1$cellID)) 
length(unique(singletCells_1$cellID))
length(unique(multiplets_1))

write.csv(singletCells_1, paste0(outputDirectory, "d2-RNA-5_Singlets.csv", sep = ""))
write.csv(multiplets_1, paste0(outputDirectory, "d2-RNA-5_Multiplets.csv"))
sample_1 <- subset(sample_1, subset = cellID %in% singletCells_1$cellID)

write.csv(singletCells_2, paste0(outputDirectory, "d5-RNA-1_Singlets.csv", sep = ""))
write.csv(multiplets_2, paste0(outputDirectory, "d5-RNA-1_Multiplets.csv"))
sample_2 <- subset(sample_2, subset = cellID %in% singletCells_2$cellID)

write.csv(singletCells_3, paste0(outputDirectory, "d5-RNA-2_Singlets.csv", sep = ""))
write.csv(multiplets_3, paste0(outputDirectory, "d5-RNA-2_Multiplets.csv"))
sample_3 <- subset(sample_3, subset = cellID %in% singletCells_3$cellID)

seed.use = 10010
#Integration     
sample_1 <- SCTransform(sample_1, verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use  = seed.use)
sample_2 <- SCTransform(sample_2, verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use  = seed.use)
sample_3 <- SCTransform(sample_3, verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use  = seed.use)


SaveSeuratRds(sample_1, file = paste0(outputRDS, "sample1_temp.rds"))
SaveSeuratRds(sample_2, file = paste0(outputRDS, "sample2_temp.rds"))
SaveSeuratRds(sample_3, file = paste0(outputRDS, "sample3_temp.rds"))
sample_1 <- LoadSeuratRds(paste0(outputRDS, "sample1_temp.rds"))
sample_2 <- LoadSeuratRds(paste0(outputRDS, "sample2_temp.rds"))
sample_3 <- LoadSeuratRds(paste0(outputRDS, "sample3_temp.rds"))

sample_1 <- RunPCA(object = sample_1, verbose = T, npcs = 100, seed.use = seed.use)
sample_2 <- RunPCA(object = sample_2, verbose = T, npcs = 100, seed.use = seed.use)
sample_3 <- RunPCA(object = sample_3, verbose = T, npcs = 100, seed.use = seed.use)

sample_list <- list(A = sample_1, 
                    B = sample_2,
                    C = sample_3)
features <- SelectIntegrationFeatures(sample_list, nfeatures = 7000)
sample_list <- PrepSCTIntegration(
  sample_list,
  anchor.features = features
)

anchors <- FindIntegrationAnchors(
  sample_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction="rpca"
)

sample <- IntegrateData(anchors, normalization.method = "SCT")

#Preprocessing
sample <- RunPCA(object = sample, verbose = T, npcs = 50, seed.use = seed.use)
sample <- FindNeighbors(object = sample, verbose = T, dims = 1:50)
sample <- FindClusters(object = sample, verbose = T, resolution = 0.5, random.seed = seed.use)
sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = 1:50)

es.max = sctype_score(scRNAseqData = sample[["SCT"]]$scale.data, scaled = TRUE, 
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

SaveSeuratRds(sample, paste0(outputRDS, "CellTagMultiHSC.rds"))
library(ggplot2)
library(ggrepel)

umap_coords <- Embeddings(sample, "umap")
df <- data.frame(umap_coords)
df$cluster <- sample$customclassif
df$cellID <- rownames(df)
write.csv(df, paste0(outputDirectoryData, "HSC.csv"))
write.csv(df, "~/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_HSC.csv")
# library(dplyr)
# library(MASS)
# 
# df$density <- with(df, MASS::kde2d(umap_1, umap_2, n = 50))
# # This function finds the highest density peak for each cell type
# find_peak <- function(data) {
#   # Compute the kernel density estimate on a grid
#   kde <- kde2d(data$umap_1, data$umap_2, n = 1000)
#   
#   # Find the grid point with the highest density
#   max_density_index <- which(kde$z == max(kde$z), arr.ind = TRUE)
#   
#   # Return the coordinates of the highest density
#   peak_coords <- cbind(kde$x[max_density_index[,1]], kde$y[max_density_index[,2]])
#   
#   # Return a data frame with peak coordinates
#   return(data.frame(umap_1 = peak_coords[,1], umap_2 = peak_coords[,2]))
# }
# 
# # Apply the function to each subset of the data frame by cell type
# peaks <- do.call(rbind, by(df, df$cluster, find_peak))
# peaks$cluster = rownames(peaks)
# plot <- DimPlot(sample, reduction = "umap", group.by = 'customclassif') + 
#   geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
#                   box.padding = unit(0.35, "lines"), 
#                   point.padding = unit(0.5, "lines"))+theme(legend.position = "none") +labs(title = NULL)
# 
# plot
# 
# ggsave(filename = paste0(outputDirectoryPlots, "HSC.png"), plot = plot, width = 6, height = 4)
# svglite(filename = paste0(outputDirectoryPlots, "HSC.svg"), width = 6, height = 4) #_withoutFiltering
# plot(plot)
# dev.off()
#  
# cellType <- as.data.frame(table(cellTypeAnnotate$Cell.Subtype))

#############################
#iEP Reprogramming
singletList <- read.csv("~/Keerthana/CellTypeCountData/AdditionalData/CellTag_FinalData/all_samples.csv")
singletList <- singletList %>% mutate(cellID = paste0(cellID, "-1"))


inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/CellTag_FinalData/matrix_files/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/CellTag/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/CellTag/"
outputRDS  <- "~/Keerthana/CellTypeCountData/finalRDS/CellTag/"
doneSample <- c()
samples <- c( "B4D12-RNA-r1-1", "B4D12-RNA-r1-2", "B4D12-RNA-r2-1", 
              "B4D12-RNA-r2-2", "B4D21-RNA-r1-1", "B4D21-RNA-r1-2", 
              "B4D21-RNA-r1-3", "B4D21-RNA-r1-4", "B4D21-RNA-r2-1", 
              "B4D21-RNA-r2-2", "B4D21-RNA-r2-3", "B4D21-RNA-r2-4", 
             "B4D3-RNA-r1-1", "B4D3-RNA-r1-2", "B4D3-RNA-r1-3", "B4D3-RNA-r1-4",
             "B4D3-RNA-r2-1", "B4D3-RNA-r2-2", "B4D3-RNA-r2-3", 
             "B4D3-RNA-r2-4", "B4D3-RNA-r2-5")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db = paste0("~/Keerthana/SingletCode_FinalData/CellTag_FinalData/scType-CellTag.xlsx")
tissue = "Embryo"
library(Matrix)
library(ggplot2)
gs_list = gene_sets_prepare(db, tissue)

seed.use = 10010
index = 0
for (i in samples){
  index = index + 1
  print(index)
  countsDir <- paste0(inputDirectory, i, "/")
  mat <- readMM(paste0(countsDir, "matrix.mtx.gz"))
  cells <- read.csv(paste0(countsDir, "barcodes.tsv"), sep = "\t", header = F)
  features <- read.csv(paste0(countsDir, "genes.tsv"), sep = "\t", header = F)
  rownames(mat) <- features$V1
  colnames(mat) <- cells$V1
  sample <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 500)
  singletCells = data.frame("cellID" = intersect(Cells(sample), singletList$cellID))
  sample$cellID = colnames(sample)
  # multiplets <- sample$cellID[!(sample$cellID %in% singletCells)]
  # write.csv(singletCells, paste0(outputDirectory, i,"_Singlets.csv", sep = ""))
  # write.csv(multiplets, paste0(outputDirectory, i, "_Multiplets.csv"))
  sample[["percent.mito"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
  sample <- filterQC(sample, 6000, 1000, 50000, 10)
  sample <- subset(sample, subset = cellID %in% singletCells$cellID)
  
  #Pre-processing
  sample <- SCTransform(sample, verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use  = seed.use)
  sample <- RunPCA(object = sample, verbose = T, npcs = 50, seed.use = seed.use)
  sample <- FindNeighbors(object = sample, verbose = T, dims = 1:50)
  sample <- FindClusters(object = sample, verbose = T, resolution = 0.5, random.seed = seed.use)
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = 1:50)
  
  es.max = sctype_score(scRNAseqData = sample[["SCT"]]$scale.data, scaled = TRUE, 
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
  
  SaveSeuratRds(sample, paste0(outputRDS,"CellTagMulti_iEP_", i,".rds"))
  library(ggplot2)
  library(ggrepel)
  
  umap_coords <- Embeddings(sample, "umap")
  df <- data.frame(umap_coords)
  df$cluster <- sample$customclassif
  df$cellID <- rownames(df) 
  write.csv(df, paste0(outputDirectoryData, i, ".csv"))
  write.csv(df, paste0("~/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_",i, ".csv" ))
  # library(dplyr)
  # library(MASS)
  # 
  # # This function finds the highest density peak for each cell type
  # find_peak <- function(data) {
  #   # Compute the kernel density estimate on a grid
  #   kde <- MASS::kde2d(data$umap_1, data$umap_2, n = 1000)
  #   
  #   # Find the grid point with the highest density
  #   max_density_index <- which(kde$z == max(kde$z), arr.ind = TRUE)
  #   
  #   # Return the coordinates of the highest density
  #   peak_coords <- cbind(kde$x[max_density_index[,1]], kde$y[max_density_index[,2]])
  #   
  #   # Return a data frame with peak coordinates
  #   return(data.frame(umap_1 = peak_coords[1, 1], umap_2 = peak_coords[1, 2]))
  # }
  # 
  # # Apply the function to each subset of the data frame by cell type
  # peaks <- do.call(rbind, by(df, df$cluster, find_peak))
  # peaks$cluster = rownames(peaks)
  # plot <- DimPlot(sample, reduction = "umap", group.by = 'customclassif') + 
  #   geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
  #                   box.padding = unit(0.35, "lines"), 
  #                   point.padding = unit(0.5, "lines"))+theme(legend.position = "none") +labs(title = NULL)
  # 
  # plot
  # 
  # ggsave(filename = paste0(outputDirectoryPlots, i, ".png"), plot = plot, width = 6, height = 4)
  # svglite(filename = paste0(outputDirectoryPlots, i, ".svg"), width = 6, height = 4) #_withoutFiltering
  # plot(plot)
  # dev.off()
  # 
  # cellType <- as.data.frame(table(df$cluster))
  # write.csv(cellType, paste0(outputDirectoryData,i,"_CellTypeCounts.csv") )
}

#Function to perform QC on a matrix and return the number of cells left and the object after QC (in case I want to save it)
filterQC <- function(sample, featuresUpper, featuresLower, countsUpper, percentMito, save = FALSE, outputPath = NULL){
  sample <- subset(x = sample, subset = nFeature_RNA > featuresLower & nFeature_RNA < featuresUpper & nCount_RNA < countsUpper & percent.mito < percentMito)
  nCells <- length(sample$cellID)
  if(save){
    SaveSeuratRds(sample, outputPath)
  }
  return(sample)
}
