lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(svglite)

cellCycleGenes <- c("Ube2c", "Hmgb2", "Hmgn2", "Tuba1b")
cellCycleGenes_2 <- c("Ccnb1", "Tubb5","Top2a", "Tubb4b")



sampleListState2 <- c("LK1_d2_1")

sampleList <- c("LK1_d2_1","LK1_d2_2", "LK1_d2_3", "LK1_d4_L_1", "LK1_d4_L_2", 
                "LK1_d4_nBC", "LK1_d4_R_1", "LK1_d4_R_4", "LK1_d6_L_1", "LK1_d6_L_2",
                "LK1_d6_R_1", "LK1_d6_R_2", "LK1_d6_R_3", "LK2_d2", "LK2_d6_1_1", "LK2_d6_2_2", "LSK_d2_3", 
                "LSK_d4_1_3", "LSK_d4_2_3", "LSK_d6_1_3", "LSK_d6_2_3",
                "LK2_d4_1", "LK2_d6_1_2", "LSK_d2_1", "LSK_d4_1_1", 
                "LSK_d4_2_1", "LSK_d6_1_1", "LSK_d6_2_1", "LK2_d4_2",
                "LK2_d6_2_1", "LSK_d2_2", "LSK_d4_1_2", "LSK_d4_2_2", 
                "LSK_d6_1_2", "LSK_d6_2_2")

singletList <- read.table(file = "~/Keerthana/CellTypeCountData/AdditionalData/LARRY_FinalData/stateFate_inVitro_metadata.txt", sep = "\t", header = T)

inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/LARRY_FinalData/"
outputDirectory <- "~/Keerthana/CellTypeCountData/LARRY/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/LARRY/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/LARRY/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/LARRY/"


######################################
#Cell Type Annotation for singlets
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db = paste0(inputDirectory, "scTypeMarker-Larry.xlsx")
tissue = "Blood"
gs_list = gene_sets_prepare(db, tissue)


library(Matrix)
mat <- readMM("~/Keerthana/CellTypeCountData/AdditionalData/LARRY_FinalData/10X/LK1_d2_1/matrix.mtx.gz")
dim(mat)

for (i in sampleList){
  countDirectory <-  paste0(inputDirectory, "10X/", i, "/")

  mat <- readMM(paste0(countDirectory, "matrix.mtx.gz"))
  mat = t(mat)
  cells <- read.csv(paste0(countDirectory, "barcodes.tsv"), sep = "\t", header = F)
  features <- read.csv(paste0(countDirectory, "genes.tsv"), sep = "\t", header = F)
  rownames(mat) <- features$V1
  colnames(mat) <- cells$V1
  sample <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200)
  #Doublet filtering
  sample$cellID = colnames(sample)
  singlets <- intersect(sample$cellID, singletList$Cell.barcode)
  multiplets <- sample$cellID[!(sample$cellID %in% singlets)]
  write.csv(singlets, paste0(outputDirectory, i, "_Singlets.csv", sep = ""))
  write.csv(multiplets, paste0(outputDirectory, i, "_Multiplets.csv"))
  sample <- subset(sample, subset = cellID %in% singletList$Cell.barcode)
  #################
  sample <- NormalizeData(sample)
  sample <- CellCycleScoring(sample, s.features = cellCycleGenes, g2m.features = cellCycleGenes_2)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 7000)
  sample <- ScaleData(sample, vars.to.regress = c("S.Score", "G2M.Score"))
  sample <- RunPCA(sample, npcs = 60)
  sample <- FindNeighbors(sample, dims = 1:60, k= 4)
  sample <- FindClusters(sample, resolution = 0.8, n.start = 500, random.seed = 1024, algorithm = 1)
  sample <- RunUMAP(sample, dims = 1:60)
  
  #Cell Type Annotation
  es.max = sctype_score(scRNAseqData = sample[["RNA"]]$scale.data, scaled = TRUE, 
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
  SaveSeuratRds(sample, paste0(outputRDS, i, ".rds"))
  
  umap_coords <- Embeddings(sample, "umap")
  df <- data.frame(umap_coords)
  df$cluster <- sample$customclassif
  colnames(df) <- c("umap_1", "umap_2", "cluster")
  write.csv(df, paste0(outputDirectoryData, i, ".csv"))
  library(dplyr)
  library(MASS)
  library(ggrepel)
  library(shadowtext)
  # This function finds the highest density peak for each cell type
  find_peak <- function(data) {
    # Compute the kernel density estimate on a grid
    kde <- MASS::kde2d(data$umap_1, data$umap_2, n = 500)
    
    # Find the grid point with the highest density
    max_density_index <- which(kde$z == max(kde$z), arr.ind = TRUE)
    
    # Return the coordinates of the highest density
    peak_coords <- cbind(kde$x[max_density_index[,1]], kde$y[max_density_index[,2]])
    
    # Return a data frame with peak coordinates
    return(data.frame(umap_1 = peak_coords[1, 1], umap_2 = peak_coords[1, 2]))
  }
  
  # Apply the function to each subset of the data frame by cell type
  peaks <- do.call(rbind, by(df, df$cluster, find_peak))
  peaks$cluster = rownames(peaks)
  plot <- DimPlot(sample, reduction = "umap", group.by = 'customclassif') + 
    geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
                    box.padding = unit(0.35, "lines"), fontface = "bold",
                    color = "black", bg.color = "white", bg.r = 0.1,
                    point.padding = unit(0.5, "lines")) +
    theme(legend.position = "none") +labs(title = NULL)
  
  plot
  
  ggsave(filename = paste0(outputDirectoryPlots, i, ".png"), plot = plot, width = 6, height = 4)
  svglite(filename = paste0(outputDirectoryPlots, i, ".svg"), width = 6, height = 4) #_withoutFiltering
  plot(plot)
  dev.off()
  
  cellType <- as.data.frame(table(df$cluster))
  write.csv(cellType, paste0(outputDirectoryData,i,"_CellTypeCounts.csv") )
}

