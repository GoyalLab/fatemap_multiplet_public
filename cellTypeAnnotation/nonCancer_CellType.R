lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(svglite)
library(ggplot2)
library(tidyverse)
#Path to directories
inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/FateMap_Heart/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/FateMapHeart/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/FateMapHeart/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/FateMapHeart/"

barcodesAll <- as.tibble(read.table(paste0("~/Keerthana/CellCounts/UnprocessedData/FateMap/non-cancer/stepFourStarcodeShavedReads50.txt"), header = TRUE))
colnames(barcodesAll) <- c("cellID", "UMI", "barcode", "sample")

sampleName_1 = "10kbarsiblingA"
sampleName_2 = "10kbarsiblingB"
seed.use = 10101
#Preprocessing the data
# Initialize the Seurat object with the raw (non-normalized data).
counts <- Read10X(data.dir = paste0(inputDirectory, "10X/", sampleName_1,"/"))
sample_1 <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
counts_2 <- Read10X(data.dir = paste0(inputDirectory, "10X/", sampleName_2,"/"))
sample_2 <- CreateSeuratObject(counts = counts_2, min.cells = 3, min.features = 200)
############
#Doublet Filtering
singletList <- doubletFiltering(barcodesAll)

singletList <- singletList %>% mutate(cellID = paste0(cellID, "-1"))
singletCells_1 = data.frame("cellID" = intersect(Cells(sample_1), singletList$cellID))
singletCells_2 = data.frame("cellID" = intersect(Cells(sample_2), singletList$cellID))

write.csv(singletCells_1, paste0(outputDirectoryData,sampleName_1,"SingletCells.csv"))
write.csv(singletCells_2, paste0(outputDirectoryData,sampleName_2,"SingletCells.csv"))


sample_1[["CellID"]] <- colnames(sample_1)
sample_2[["CellID"]] <- colnames(sample_2)
#get the required cells
#subset the required cells
sample_1 <- subset(sample_1, subset = CellID %in% singletCells_1$cellID)
sample_2 <- subset(sample_2, subset = CellID %in% singletCells_2$cellID)

####################
#Sample QC


#subset based on feature cutoffs
featuresUpper = 10000
sample_1 <- subset(sample_1, cells = colnames(sample_1)[which(sample_1$nFeature_RNA < featuresUpper)])
sample_2 <- subset(sample_2, cells = colnames(sample_2)[which(sample_2$nFeature_RNA < featuresUpper)])
message("Filtered dimensions")
print(dim(sample_1))
print(dim(sample_2))

sample_1 <- SCTransform(sample_1, verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use  = seed.use)
sample_2 <- SCTransform(sample_2, verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use  = seed.use)

sample_list <- list(A = sample_1, 
                    B = sample_2)
features <- SelectIntegrationFeatures(sample_list, nfeatures = 5000)
sample_list <- PrepSCTIntegration(
  sample_list,
  anchor.features = features
)

anchors <- FindIntegrationAnchors(
  sample_list,
  normalization.method = "SCT",
  anchor.features = features
)

sample <- IntegrateData(anchors, normalization.method = "SCT")

#Dimensional reduction & clustering
message("Running dimensional reduction algos")

sample <- RunPCA(object = sample, verbose = T, npcs = 50, seed.use = seed.use)

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
db = paste0(inputDirectory, "/FateMapHeart-scType.xlsx")
tissue = "Heart"
gs_list = gene_sets_prepare(db, tissue)

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

SaveSeuratRds(sample, paste0(outputRDS, "FateMapHeart", ".rds"))

umap_coords <- Embeddings(sample, "umap")
df <- data.frame(umap_coords)
df$cluster <- sample$customclassif
colnames(df) <- c("umap_1", "umap_2", "cluster")
df$cellID <- rownames(df)
write.csv(df, paste0(outputDirectoryData, "FateMapHeart", ".csv"))
write.csv(df, "~/Keerthana/CellTypeCountData/plotData/UmapData/Jiang_et_al/Jiang_et_al.csv")
# library(dplyr)
# library(MASS)
# library(ggrepel)
# library(shadowtext)
# # This function finds the highest density peak for each cell type
# find_peak <- function(data) {
#   # Compute the kernel density estimate on a grid
#   kde <- MASS::kde2d(data$umap_1, data$umap_2, n = 500)
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
#                   box.padding = unit(0.35, "lines"), fontface = "bold",
#                   color = "black", bg.color = "white", bg.r = 0.1,
#                   point.padding = unit(0.5, "lines")) +
#   theme(legend.position = "none") +labs(title = NULL)
# 
# plot
# 
# ggsave(filename = paste0(outputDirectoryPlots, "FateMapHeart", ".png"), plot = plot, width = 6, height = 4)
# svglite(filename = paste0(outputDirectoryPlots, "FateMapHeart", ".svg"), width = 6, height = 4) #_withoutFiltering
# plot(plot)
# dev.off()
# 
# cellType <- as.data.frame(table(df$cluster))
# write.csv(cellType, paste0(outputDirectoryData,"FateMapHeart","_CellTypeCounts.csv") )

#Singlet filtering method from the paper
doubletFiltering <- function(barcodes10X){
  umiCut = 2 # Minimum UMI cutoff (per barcode) for reliable analysis. NOTE that unless this is taking place after unique() subset of barcodes it is not actually an umiCut but rather cutoff for total number of reads (without accounting for PCR duplicates)
  fracumiCut = 0.3
  #to a point, increasing fracumiCut decreases total cells recovered but increases proportion with 1 barcode
  
  linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.
  
  barcodePostUmiCutoff = barcodes10X %>% 
    group_by(cellID, barcode, sample) %>%
    summarise(nUMI = length(sample)) %>%
    filter(nUMI >= umiCut) %>% 
    group_by(cellID,sample) %>%
    mutate(fracUMI = nUMI/sum(nUMI)) %>%
    filter(fracUMI >= fracumiCut) %>%
    group_by(cellID, sample) %>%
    mutate(nLineages = length(cellID)) 
  #6952 for umiCut=2+fracumiCut=0.3
  
  ###taking only single cellid-->barcode mappings; removes doublets, ambient RNA (more than 1 barcode per cellID)
  barcodePostUmiLineageCutoff = barcodePostUmiCutoff %>%
    ungroup() %>%
    filter(nLineages <= linCut) %>%
    unique() 
  return(barcodePostUmiLineageCutoff)
}

