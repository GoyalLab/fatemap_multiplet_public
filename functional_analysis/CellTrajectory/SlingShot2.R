library(Seurat)
library(readr)
library(dplyr)
library(slingshot)
set.seed(10101)
library(svglite)
library(ggplot2)
library(DelayedMatrixStats)
library(igraph)
library(tradeSeq)
#Need matrix version 1.1
library(matrixStats)
inputDirectoryAll <- "~/Keerthana/singletCode/CellTrajectory/"
plotPNG <- "~/Keerthana/CellTrajectory/Plots/"
alreadyRun <- c("FM01", "FM02", "FM03", "FM04", "FM05",
                "FM06", "FM08",  "Biorxiv", "cellTag",
                "ClonMapper", "LARRY", "non_cancer", "smartseq3_umis", "watermelon")
not_run <- c("smartseq3_reads", "TREX")
to_run <- c("FM04")
experimentID <- c("TREX")
samples <- c("control","ten", "twenty", "forty")

for (experimentID in to_run){
  inputDirectory <- paste0(inputDirectoryAll, experimentID, "/data/")
  print(experimentID)
  for (sampleName in samples) {
    dataPath <- paste0(inputDirectory, sampleName, ".rds")
    print(sampleName)
    sample <- readRDS(dataPath)
    sample$clusterMod <- droplevels(sample$seurat_clusters)
    print("CP1")
    reducedDimMatrix <- sample@reductions$pca@cell.embeddings[,1:10 ]
    clusterList <- sample$clusterMod
    PCA <- as.data.frame(reducedDimMatrix)
    PCA$cluster <- as.vector(clusterList)
    PCA$singlet <- sample$label
    rm(sample) # Consider removing if you need 'sample' for later use
    gc()
    print("CP2")
    lineage <- getLineages(reducedDimMatrix, clusterList)
    print("CP3")
    curves <- getCurves(lineage)
    curvesPlot <- slingCurves(curves, as.df = TRUE)
    # write_csv(curvesPlot,file = paste0(inputDirectoryAll, "data/curvesPlot/", experimentID, "_", sampleName, ".csv") )
    PCA$cellID <- rownames(PCA)
    pseudoTime <- data.frame(curves@assays@data@listData$pseudotime)
    # write_csv(pseudoTime, file = paste0(inputDirectoryAll, "data/pseudotime/", experimentID, "_", sampleName, ".csv"))
    cellWeights <- data.frame(curves@assays@data@listData$weights)
    # write_csv(cellWeights, file = paste0(inputDirectoryAll, "data/cellWeights/", experimentID, "_", sampleName, ".csv"))
    
    mst <- slingMST(lineage)
    # write_graph(mst, file = paste0(inputDirectoryAll, "lineageGraph/", experimentID, "_", sampleName, ".csv"))
    # write_csv(PCA, file = paste0(inputDirectoryAll, "data/PCA/", experimentID, "_", sampleName, ".csv"))
    
    lineages_txt <- sapply(lineage@metadata$lineages, function(x) paste(x, collapse = ", "))
    lineages_df <- data.frame(Lineage = lineages_txt, stringsAsFactors = FALSE)
    # write.csv(lineages_df, paste0(inputDirectoryAll, "lineageTrajectory/", experimentID, "_", sampleName, ".csv"), row.names = FALSE, quote = FALSE) # Fixed parentheses and arguments
    
    plotDoublets <- ggplot(PCA, aes(x = PC_1, y = PC_2)) + # Assuming PC1, PC2 are correct column names
      geom_point(data = PCA[PCA$singlet == "singlet",], aes(x = PC_1, y = PC_2), color = "grey", alpha = 0.8) +
      geom_point(data = PCA[PCA$singlet == "doublet",], aes(x = PC_1, y = PC_2), color = "red3", alpha = 0.6) +
      geom_path(data = curvesPlot %>% arrange(Order), aes(group = Lineage)) +
      theme_void() +
      theme(plot.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(plotDoublets)
    plotClusters <- ggplot(PCA, aes(x = PC_1, y = PC_2)) + 
      geom_point(aes(color = factor(cluster)), alpha = 0.6) + 
      geom_path(data= curvesPlot %>% arrange(Order), aes(group = Lineage)) + 
      theme_void() +
      theme(legend.position = "none",
            plot.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    plotClusters
    # svglite(file = paste0(plotPNG, "doubletsSVG/", experimentID, "_", sampleName, "_doublets.svg"), width = 4, height = 6)
    # print(plotDoublets) # Use print to ensure plots are rendered in non-interactive sessions
    # dev.off()
    
    pngFilePathCluster <- paste0(plotPNG, "clusterPNG/", experimentID, "_", sampleName, "_cluster.png")
    # ggsave(plotClusters, file = pngFilePathCluster, width = 6, height = 4)
    
    pngFilePathDoublets <- paste0(plotPNG, "doubletsPNG/", experimentID, "_", sampleName, "_doublets.png")
    # ggsave(plotDoublets, file = pngFilePathDoublets, width = 4, height = 6)
  }
}


############################
#Trying out tradeSeq
set.seed(5)
# Assuming 'sample' is your merged Seurat object


library(BiocParallel)
bpparam <- SnowParam(type = "SOCK", progressbar = TRUE)

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 40

counts <- sample@assays$RNA$counts[, 1:1000]
icMat <- evaluateK(counts = raw_counts, sds = curves, k = 5:10, nGenes = 200, verbose = T, parallel = T, BPPARAM = BPPARAM)
experimentID = "smartseq3_umis"
experimentTraj <- c("FM06")
for (experimentID in experimentTraj){
  inputDirectory <- paste0(inputDirectoryAll, experimentID, "/data/")
  print(experimentID)
  for (sampleName in samples){
    dataPath <- paste0(inputDirectory, sampleName, ".rds")
    print(sampleName)
    sample <- readRDS(dataPath)
    reducedDimMatrix <- sample@reductions$pca@cell.embeddings[, 1:10]
    clusterList <- sample$seurat_clusters
    PCA <- as.data.frame(reducedDimMatrix)
    PCA$cluster <- as.vector(clusterList)
    PCA$singlet <- sample$label
    lineage <- getLineages(reducedDimMatrix, clusterList)
    curves <- getCurves(lineage)
    
    sample <- JoinLayers(sample)
    variableGenes <- VariableFeatures(sample)
    
    raw_counts <- GetAssayData(object = sample, assay = "RNA", layer = "counts")[variableGenes,]
    sce <- fitGAM(counts = as.matrix(raw_counts), sds = curves, parallel = F, BPPARAM = BPPARAM)
    saveRDS(sce, paste0("~/Keerthana/CellTrajectory/sceTrajectory/", experimentID,"_", sampleName, ".rds"))
    endRes <- diffEndTest(sce)
    endRes$Gene <- rownames(endRes)
    startEndTest <- startVsEndTest(sce)
    startEndTest$Gene <- rownames(startEndTest)
    write_csv(startEndTest, paste0("~/Keerthana/CellTrajectory/sceTrajectory/sVE_", experimentID,"_", sampleName, ".csv"))
    write_csv(endRes, paste0("~/Keerthana/CellTrajectory/sceTrajectory/DEG_", experimentID,"_", sampleName, ".csv"))
    }
}


