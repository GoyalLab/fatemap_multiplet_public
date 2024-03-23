lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(ggplot2)
library(svglite)
#Path to directories

library(reticulate)

if (!(py_module_available("scanorama"))){
  message("Installing 'Scanorama' python module...")
  py_install("scanorama", pip = T)
}
scanorama <- import('scanorama')

inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/FateMap_FinalData/NJ01/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/FateMap-NonCancer/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/FateMap-NonCancer/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/FateMap-NonCancer/"

barcodesAll <- as.tibble(read.table("~/Keerthana/CellCounts/UnprocessedData/FateMap/NJ01/stepFourStarcodeShavedReads50.txt", header = TRUE))
colnames(barcodesAll) <- c("cellID", "barcode", "sample")
samples = c("filtered_feature_bc_matrix", "filtered_feature_bc_matrix2", 
            "filtered_feature_bc_matrix3", "filtered_feature_bc_matrix4",
            "filtered_feature_bc_matrix5", "filtered_feature_bc_matrix6")
group_names <- c("DMSO_A", "DMSO_B", "LSD1i_A", "LSD1i_B", "DOT1Li_A", "DOT1Li_B")
s1_data <- Read10X(data.dir = paste0(inputDirectory, "10X/filtered_feature_bc_matrix"))
s2_data <- Read10X(data.dir = paste0(inputDirectory, "10X/filtered_feature_bc_matrix2"))
s3_data <- Read10X(data.dir = paste0(inputDirectory, "10X/filtered_feature_bc_matrix3"))
s4_data <- Read10X(data.dir = paste0(inputDirectory, "10X/filtered_feature_bc_matrix4"))
s5_data <- Read10X(data.dir = paste0(inputDirectory, "10X/filtered_feature_bc_matrix5"))
s6_data <- Read10X(data.dir = paste0(inputDirectory, "10X/filtered_feature_bc_matrix6"))

s1 <- CreateSeuratObject(counts = s1_data, project = "DMSO_A", min.cells = 3, min.features = 200)
s2 <- CreateSeuratObject(counts = s2_data, project = "DMSO_B", min.cells = 3, min.features = 200)
s3 <- CreateSeuratObject(counts = s3_data, project = "LSD1i_A", min.cells = 3, min.features = 200)
s4 <- CreateSeuratObject(counts = s4_data, project = "LSD1i_B", min.cells = 3, min.features = 200)
s5 <- CreateSeuratObject(counts = s5_data, project = "DOT1Li_A", min.cells = 3, min.features = 200)
s6 <- CreateSeuratObject(counts = s6_data, project = "DOT1Li_B", min.cells = 3, min.features = 200)

#Doublet Filtering
#Doublet Filtering
singletList <- doubletFiltering(barcodesAll)
singletList$cellID <- paste0(singletList$cellID, "-1")
singletCells_1 = data.frame("cellID" = intersect(Cells(s1), singletList$cellID))
singletCells_2 = data.frame("cellID" = intersect(Cells(s2), singletList$cellID))
singletCells_3 = data.frame("cellID" = intersect(Cells(s3), singletList$cellID))
singletCells_4 = data.frame("cellID" = intersect(Cells(s4), singletList$cellID))
singletCells_5 = data.frame("cellID" = intersect(Cells(s5), singletList$cellID))
singletCells_6 = data.frame("cellID" = intersect(Cells(s6), singletList$cellID))

write.csv(singletCells_1, paste0(outputDirectoryData,group_names[1],"_SingletCells.csv"))
write.csv(singletCells_2, paste0(outputDirectoryData,group_names[2],"_SingletCells.csv"))
write.csv(singletCells_3, paste0(outputDirectoryData,group_names[3],"_SingletCells.csv"))
write.csv(singletCells_4, paste0(outputDirectoryData,group_names[4],"_SingletCells.csv"))
write.csv(singletCells_5, paste0(outputDirectoryData,group_names[5],"_SingletCells.csv"))
write.csv(singletCells_6, paste0(outputDirectoryData,group_names[6],"_SingletCells.csv"))
# 
# write.csv(multipletCells_1, paste0(outputDirectory,group_names[1],"MultipletCells.csv"))
# write.csv(multipletCells_2, paste0(outputDirectory,group_names[2],"MultipletCells.csv"))
# write.csv(multipletCells_3, paste0(outputDirectory,group_names[3],"MultipletCells.csv"))
# write.csv(multipletCells_4, paste0(outputDirectory,group_names[4],"MultipletCells.csv"))
# write.csv(multipletCells_5, paste0(outputDirectory,group_names[5],"MultipletCells.csv"))
# write.csv(multipletCells_6, paste0(outputDirectory,group_names[6],"MultipletCells.csv"))

s1[["CellID"]] <- colnames(s1)
s2[["CellID"]] <- colnames(s2)
s3[["CellID"]] <- colnames(s3)
s4[["CellID"]] <- colnames(s4)
s5[["CellID"]] <- colnames(s5)
s6[["CellID"]] <- colnames(s6)
#get the required cells
#subset the required cells
s1 <- subset(s1, subset = CellID %in% singletCells_1$cellID)
s2 <- subset(s2, subset = CellID %in% singletCells_2$cellID)
s3 <- subset(s3, subset = CellID %in% singletCells_3$cellID)
s4 <- subset(s4, subset = CellID %in% singletCells_4$cellID)
s5 <- subset(s5, subset = CellID %in% singletCells_5$cellID)
s6 <- subset(s6, subset = CellID %in% singletCells_6$cellID)

#Preprocessing-data
s1[["percent.mt"]] <- PercentageFeatureSet(object = s1, pattern = "^MT-")
s2[["percent.mt"]] <- PercentageFeatureSet(object = s2, pattern = "^MT-")
s3[["percent.mt"]] <- PercentageFeatureSet(object = s3, pattern = "^MT-")
s4[["percent.mt"]] <- PercentageFeatureSet(object = s4, pattern = "^MT-")
s5[["percent.mt"]] <- PercentageFeatureSet(object = s5, pattern = "^MT-")
s6[["percent.mt"]] <- PercentageFeatureSet(object = s6, pattern = "^MT-")

s1 <- subset(x = s1, subset = nFeature_RNA > 200 & percent.mt < 30)
s2 <- subset(x = s2, subset = nFeature_RNA > 200 & percent.mt < 30)
s3 <- subset(x = s3, subset = nFeature_RNA > 200 & percent.mt < 30)
s4 <- subset(x = s4, subset = nFeature_RNA > 200 & percent.mt < 30)
s5 <- subset(x = s5, subset = nFeature_RNA > 200 & percent.mt < 30)
s6 <- subset(x = s6, subset = nFeature_RNA > 200 & percent.mt < 30)

scTransform <- merge(s1, y = c(s2, s3, s4, s5, s6), add.cell.ids = c("S1", "S2", "S3", "S4", "S5", "S6"))



# Repeat the logic for different group_names
group_names <- c("DMSO_A", "DMSO_B", "LSD1i_A", "LSD1i_B", "DOT1Li_A", "DOT1Li_B")

datasets <- list()
gene_lists <- list()

for (group_name in group_names) {
  # Get the expression data for the specified assay
  assay_data <- scTransform[["RNA"]][group_name]
  print(dim(assay_data))
  # Transpose the data matrix (convert to matrix if not already)
  transposed_data <- t(as.matrix(assay_data))
  print(dim(transposed_data))
  group_indices <- which(Idents(scTransform) == group_name)
  # Store results in lists
  datasets[[group_name]] <- transposed_data
  gene_lists[[group_name]] <- rownames(assay_data)
}

names(datasets) = NULL
names(gene_lists) = NULL

integrated_corrected_data = scanorama$correct(datasets, gene_lists, return_dimred = TRUE, return_dense = TRUE, ds_names = group_names, verbose = TRUE)

corrected_scanorama <- t(do.call(rbind, integrated_corrected_data[[2]]))
colnames(corrected_scanorama) <- colnames(scTransform)
rownames(corrected_scanorama) <- integrated_corrected_data[[3]]
corrected_scanorama_pca <- t(do.call(rbind, integrated_corrected_data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(scTransform)

scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
scTransform[["scanorama"]] <- scanorama_assay
DefaultAssay(scTransform) <- "scanorama"



all.genes <- rownames(scTransform)
scTransform <- ScaleData(scTransform, features = all.genes)
scTransform <- RunPCA(object = scTransform, assay = "scanorama", features= all.genes, reduction.name = "pca_scanorama")

scTransform <- FindNeighbors(object=scTransform, features = all.genes, dims = 1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
scTransform <- FindClusters(object=scTransform,graph.name = "scanorama_snn", resolution = 0.5)
scTransform <- RunUMAP(object = scTransform, reduction = "pca_scanorama", dims = 1:50, reduction.name = "umap_scanorama")


######################################
#Cell Type Annotation for singlets
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db = paste0(inputDirectory, "/scTypeMarker_fateMapNonCancer.xlsx")
tissue = "stemCell"
gs_list = gene_sets_prepare(db, tissue)

es.max = sctype_score(scRNAseqData = scTransform[["scanorama"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = NULL) 

cL_resutls = do.call("rbind", lapply(unique(scTransform@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(scTransform@meta.data[scTransform@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scTransform@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

scTransform@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  scTransform@meta.data$customclassif[scTransform@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

SaveSeuratRds(scTransform, paste0(outputRDS, "FateMap-NonCancer", ".rds"))

umap_coords <- Embeddings(scTransform, "umap_scanorama")
df <- data.frame(umap_coords)
df$cluster <- scTransform$customclassif
colnames(df) <- c("umap_1", "umap_2", "cluster")
df$cellID <- rownames(df)
write.csv(df, paste0(outputDirectoryData, "FateMap-NonCancer", ".csv"))
write_csv(df, paste0("~/Keerthana/CellTypeCountData/plotData/UmapData/Jain_et_al/jain_et_al.csv"))
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
plot <- DimPlot(scTransform, reduction = "umap_scanorama", group.by = 'customclassif') + 
  geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
                  box.padding = unit(0.35, "lines"), fontface = "bold",
                  color = "black", bg.color = "white", bg.r = 0.1,
                  point.padding = unit(0.5, "lines")) +
  theme(legend.position = "none") +labs(title = NULL)

plot

ggsave(filename = paste0(outputDirectoryPlots, "FateMap-NonCancer", ".png"), plot = plot, width = 6, height = 4)
svglite(filename = paste0(outputDirectoryPlots, "FateMap-NonCancer", ".svg"), width = 6, height = 4) #_withoutFiltering
plot(plot)
dev.off()

cellType <- as.data.frame(table(df$cluster))
write.csv(cellType, paste0(outputDirectoryData,"FateMap-NonCancer","_CellTypeCounts.csv") )

#Singlet filtering method from the paper
doubletFiltering <- function(barcodes10X){
  
  #Checking for pattern for 10X barcodes - ideally all of the barcodes in this file must be 10X - so a check
  barcodes10X %>% .$cellID %>% unique() %>% length()
  barcodes10XFilter <- barcodes10X[grepl("([AT][CG][ATCG]){8,}", barcodes10X$barcode), ]
  barcodes10XFilter %>% .$cellID %>% unique() %>% length()
  
  #Filtering barcodes based on UMI counts and percentage of UMI counts from a barcode vs total (< 0.4 means it is removed)
  umiCut = 0
  fracUMICut = 0.4
  
  upToLineageCounting <- barcodes10XFilter %>% group_by(cellID, barcode, sample) %>% summarise(nUMI = length(cellID)) %>% filter(nUMI >= umiCut) %>%
    group_by(cellID, sample) %>% mutate(fracUMI = nUMI/sum(nUMI)) %>% filter(fracUMI >= fracUMICut) %>%
    group_by(cellID) %>% mutate(nLineages = length(cellID))
  
  upToLineageCounting %>% .$cellID %>% unique() %>% length()
  
  #Removing cells containing more than one barcode even after the above filtering
  linCut = 1
  linCountToOverlaps <- upToLineageCounting %>% ungroup() %>% filter(nLineages <= linCut) %>% unique()
  linCountToOverlaps %>% .$cellID %>% length()
  return(linCountToOverlaps)
}
