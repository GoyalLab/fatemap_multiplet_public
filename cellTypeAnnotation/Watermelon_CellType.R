lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(svglite)
#Path to directories
inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/Watermelon10X_FinalData/"
outputDirectory <- "~/Keerthana/CellTypeCountData/Watermelon/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/Watermelon/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/Watermelon/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/Watermelon/"
sampleName = "T47D-naive-2"
seed.use = 10101
#Preprocessing the data

sample <- readRDS("~/Keerthana/RawData/Watermelon_RawData/T47D_seurat_scRNAseq.rds")
# Initialize the Seurat object with the raw (non-normalized data).
counts <- Read10X(data.dir = paste0(inputDirectory, sampleName,"/"))
sample <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

#Sample QC

# calculate mito percentage
sample[["percent.mito"]] <- PercentageFeatureSet(sample, pattern = "^MT-")

#subset based on feature cutoffs and mito genes percentage
featuresUpper = 8000
featuresLower = 2500
mitoPercentage = 25
countsUpper = 60000
summary(sample$percent.mito)
summary(sample$nFeature_RNA)
sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA < countsUpper)])
sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA > featuresLower)])
sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA < featuresUpper)])
sample <- subset(sample, cells = colnames(sample)[which(sample$percent.mito < mitoPercentage)])

sample <- NormalizeData(sample)
#Regressing out cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sample$CC.Difference <- sample$S.Score - sample$G2M.Score

sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 5000)
sample <- ScaleData(sample, vars.to.regress = c("CC.Difference"), block.size = 3000, features = )

sample <- RunPCA(object = sample, verbose = T, npcs = 50, seed.use = seed.use)

# find neighbors based on 1 to 50 dimensions
sample <- FindNeighbors(object = sample, verbose = T, dims = 1:50)

# Run tSNE & UMAP on 1 to dims
sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = 1:50)
sample <- FindClusters(object = sample, verbose = T, resolution = 0.8, random.seed = seed.use)


######################################
#Cell Type Annotation for singlets
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db = paste0(inputDirectory, "/Watermelon10X_scTypeMarkers.xlsx")
tissue = "Breast"
gs_list = gene_sets_prepare(db, tissue)

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
SaveSeuratRds(sample, paste0(outputRDS, "Watermelon.rds"))

sample <- readRDS(paste0(outputRDS,"T47D.rds"))
umap_coords <- Embeddings(sample, "umap")

df <- read.csv(paste0(outputDirectoryData, "T47D.csv"))

df <- data.frame(umap_coords)
df$cluster <- sample$customclassif
colnames(df) <- c("umap_1", "umap_2", "cluster")
df$cellID <- rownames(df)
write.csv(df, paste0(outputDirectoryData, "T47D.csv"))
write_csv(df, "~/Keerthana/CellTypeCountData/plotData/UmapData/Watermelon/Watermelon_T47D.csv")
library(dplyr)
library(MASS)
library(ggrepel)
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
                  color = "black", bg.color = "white", bg.r = 0.1, nudge_x = -0.5,
                  point.padding = unit(0.5, "lines")) +
  theme(legend.position = "none") +labs(title = NULL)

plot

ggsave(filename = paste0(outputDirectoryPlots, "T47D.png"), plot = plot, width = 6, height = 4)
svglite(filename = paste0(outputDirectoryPlots, "T47D.svg"), width = 6, height = 4) #_withoutFiltering
plot(plot)
dev.off()

cellType <- as.data.frame(table(df$cluster))
write.csv(cellType, paste0(outputDirectoryData,"T47D_CellTypeCounts.csv") )

df_umap <- as.data.frame(Watermelon@reductions$umap@cell.embeddings)
df_umap$cellID <- rownames(df_umap)
df_umap$cluster <- Watermelon$customclassif
write_csv(df_umap, "~/Keerthana/CellTypeCountData/plotData/UmapData/Watermelon/Watermelon_T47D.csv")
