lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(svglite)
library(reticulate)
if (!(py_module_available("scanorama"))){
  message("Installing 'Scanorama' python module...")
  py_install("scanorama", pip = T)
}
scanorama <- import('scanorama')

#Path to directories
inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/FateMap_FinalData/"
outputDirectory <- "~/Keerthana/CellTypeCountData/FateMap/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/FateMap/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/FateMap/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/FateMap/"
singletListAll <- read.csv("~/Keerthana/CellTypeCountData/AdditionalData/FateMap_FinalData/singletList.csv")
experimentID = "FM01"

#Used for FM01
s1_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/sample1"))
s2_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/sample2"))
s3_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/sample3"))
s4_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/sample4"))

#used for FM02, {FM03, FM04, FM05}(only 1 and 2)
s1_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/filtered_feature_bc_matrix"))
s2_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/filtered_feature_bc_matrix2"))
s3_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/filtered_feature_bc_matrix3"))
s4_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/filtered_feature_bc_matrix4"))

s1_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/FM06-WM989Naive-1/filtered_feature_bc_matrix"))
s2_data <- Read10X(data.dir=paste0(inputDirectory, experimentID,"/10X/FM06-WM989Naive-2/filtered_feature_bc_matrix"))

s1 <- CreateSeuratObject(counts = s1_data, min.cells = 3, min.features = 200, project = "S1")
s2 <- CreateSeuratObject(counts = s2_data, min.cells = 3, min.features = 200, project = "S2")
s3 <- CreateSeuratObject(counts = s3_data, min.cells = 3, min.features = 200, project = "S3")
s4 <- CreateSeuratObject(counts = s4_data, min.cells = 3, min.features = 200, project = "S4")

#

singletCells_1 = data.frame("cellID" = intersect(Cells(s1), singletList$V1))
singletCells_2 = data.frame("cellID" = intersect(Cells(s2), singletList$V1))
singletCells_3 = data.frame("cellID" = intersect(Cells(s3), singletList$V1))
singletCells_4 = data.frame("cellID" = intersect(Cells(s4), singletList$V1))

s1[["CellID"]] <- colnames(s1)
s2[["CellID"]] <- colnames(s2)
s3[["CellID"]] <- colnames(s3)
s4[["CellID"]] <- colnames(s4)

#subset the required cells
s1 <- subset(s1, subset = CellID %in% singletCells_1$cellID)
s2 <- subset(s2, subset = CellID %in% singletCells_2$cellID)
s3 <- subset(s3, subset = CellID %in% singletCells_3$cellID)
s4 <- subset(s4, subset = CellID %in% singletCells_4$cellID)

#Preprocessing-data
s1[["percent.mt"]] <- PercentageFeatureSet(object = s1, pattern = "^MT-")
s2[["percent.mt"]] <- PercentageFeatureSet(object = s2, pattern = "^MT-")
s3[["percent.mt"]] <- PercentageFeatureSet(object = s3, pattern = "^MT-")
s4[["percent.mt"]] <- PercentageFeatureSet(object = s4, pattern = "^MT-")


s1 <- subset(x = s1, subset = nFeature_RNA < 7200 & percent.mt < 26)
s2 <- subset(x = s2, subset = nFeature_RNA < 7200 & percent.mt < 26)
s3 <- subset(x = s3, subset = nFeature_RNA < 7200 & percent.mt < 26)
s4 <- subset(x = s4, subset = nFeature_RNA < 7200 & percent.mt < 26)

#Integration using scanorama
sample_merged <- merge(s1, y = c(s2, s3, s4), add.cell.ids = c("S1", "S2", "S3", "S4"))
group_names <- c("S1", "S2", "S3", "S4")
sample_merged <- merge(s1, y = c(s2, s3), add.cell.ids = c("S1", "S2", "S3"))
group_names <- c("S1", "S2", "S3")
sample_merged <- merge(s1, y = c(s2), add.cell.ids = c("S1", "S2"))
group_names <- c("S1", "S2")
datasets <- list()
gene_lists <- list()

for (group_name in group_names) {
  # Get the expression data for the specified assay
  assay_data <- sample_merged[["RNA"]][group_name]
  print(dim(assay_data))
  # Transpose the data matrix (convert to matrix if not already)
  transposed_data <- t(as.matrix(assay_data))
  print(dim(transposed_data))
  group_indices <- which(Idents(sample_merged) == group_name)
  # Store results in lists
  datasets[[group_name]] <- transposed_data
  gene_lists[[group_name]] <- rownames(assay_data)
}

names(datasets) = NULL
names(gene_lists) = NULL

integrated_corrected_data = scanorama$correct(datasets, gene_lists, return_dimred = TRUE, return_dense = TRUE, ds_names = group_names, verbose = TRUE)

corrected_scanorama <- t(do.call(rbind, integrated_corrected_data[[2]]))
colnames(corrected_scanorama) <- colnames(sample_merged)
rownames(corrected_scanorama) <- integrated_corrected_data[[3]]
corrected_scanorama_pca <- t(do.call(rbind, integrated_corrected_data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(sample_merged)

scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
sample_merged[["scanorama"]] <- scanorama_assay
DefaultAssay(sample_merged) <- "scanorama"

#Pre-processing
all.genes <- rownames(sample_merged)
sample_merged <- ScaleData(sample_merged, features = all.genes)
sample_merged <- RunPCA(object = sample_merged, assay = "scanorama", features= all.genes, reduction.name = "pca_scanorama", verbose = TRUE)

sample_merged <- FindNeighbors(object=sample_merged, features = all.genes, dims = 1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
sample_merged <- FindClusters(object=sample_merged,graph.name = "scanorama_snn", resolution = 0.8)
sample_merged <- RunUMAP(object = sample_merged, reduction = "pca_scanorama", dims = 1:50, reduction.name = "umap_scanorama")


######################################
#Cell Type Annotation for singlets
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db = paste0("~/Keerthana/CellTypeCountData/AdditionalData/FateMap_FinalData/FateMapMarkers.xlsx")
tissue = "Melanoma"

# db = paste0("~/Keerthana/CellTypeCountData/AdditionalData/Watermelon10X_FinalData/Watermelon10X_scTypeMarkers.xlsx")
# tissue = "Breast"

gs_list = gene_sets_prepare(db, tissue)

es.max = sctype_score(scRNAseqData = sample_merged[["scanorama"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = NULL) 

cL_resutls = do.call("rbind", lapply(unique(sample_merged@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sample_merged@meta.data[sample_merged@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sample_merged@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

sample_merged@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sample_merged@meta.data$customclassif[sample_merged@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

SaveSeuratRds(sample_merged, paste0(outputRDS, experimentID, ".rds"))

umap_coords <- Embeddings(sample_merged, "umap_scanorama")
df <- data.frame(umap_coords)
df$cluster <- sample_merged$customclassif
colnames(df) <- c("umap_1", "umap_2", "cluster")
write.csv(df, paste0(outputDirectoryData, experimentID, ".csv"))
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
plot <- DimPlot(sample_merged, reduction = "umap_scanorama", group.by = 'customclassif') + 
  geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
                  box.padding = unit(0.35, "lines"), fontface = "bold",
                  color = "black", bg.color = "white", bg.r = 0.1,
                  point.padding = unit(0.5, "lines")) +
                  theme(legend.position = "none") +labs(title = NULL)

plot

ggsave(filename = paste0(outputDirectoryPlots, experimentID, ".png"), plot = plot, width = 6, height = 4)
svglite(filename = paste0(outputDirectoryPlots, experimentID, ".svg"), width = 6, height = 4) #_withoutFiltering
plot(plot)
dev.off()

cellType <- as.data.frame(table(df$cluster))
write.csv(cellType, paste0(outputDirectoryData,experimentID,"_CellTypeCounts.csv") )

