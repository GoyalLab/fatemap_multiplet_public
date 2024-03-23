# runSeuratPipeline
# Dane Vassiliadis
# 8-6-20

# Pipeline to automate basic 10X sample processing and visualisation using
# Seurat within R


#-----------------------
# Required libs
#-----------------------
# 
# remotes::install_version("matrixStats", version="1.1.0")
library(matrixStats)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(dplyr)
# BiocManager::install("PCAtools")
library(stringr)
# BiocManager::install(version = "3.16")
library(PCAtools)
library(HGNChelper)

#Cell cycle genes
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

g2m.genes.mouse <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a",
                     "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf",
                     "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb",
                     "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp",
                     "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1",
                     "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr",
                     "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf",
                     "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")
s.genes.mouse <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4",
                   "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl",
                   "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rpa2",
                   "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2",
                   "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2",
                   "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm",
                   "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1",
                   "E2f8")
#-----------------------
# Main function
#-----------------------
seed.use  = 10101
dims = 1:30
inputdir = "~/Keerthana/CellTypeCountData/AdditionalData/SPLINTR_FinalData/"
outdir = "~/Keerthana/CellTypeCountData/SPLINTR/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/SPLINTR/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/SPLINTR/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/SPLINTR/"

input.dir <- file.path(inputdir)
output.dir <- file.path(outdir)

samples <- c("retransplant","chemoDay2_1", "chemoDay2_2", "chemoDay5_1", "chemoDay5_2", "chemoDay7_1",
             "chemoDay7_2", "chemoKRASBaseline_1", "chemoKRASBaseline_2", "chemoVehicle_1",
             "chemoVehicle_2", "inVitro_FLT3", "inVitro_KRAS", "inVitro_MLL_1", "inVitro_MLL_2")


sampleName <- c("inVitro_MLL_2")
for (i in samples){
  Processing_func(i)
}
Processing_func <- function(sampleName){
  counts.dir = paste0("~/Keerthana/SingletCode_FinalData/SPLINTR_FinalData/10X/", sampleName, "/")
  counts <- Read10X(counts.dir, strip.suffix = T)
  sample <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
  
  ############
  #Doublet Filtering
  
  barcodeCount <- read.csv(paste0(input.dir,"all_samples_barcode_counts.csv"))
  barcodeCount = barcodeCount[barcodeCount$sample == sampleName,]
  result <- barcodeCount %>%
    group_by(cellID, barcode) %>%
    summarize(n = n())
  
  fin <- result %>%
    group_by(cellID) %>%
    summarize(n = n())
  fin <- fin %>%
    mutate(doubletState = ifelse(n == 1, "singlet", "doublet"))
  # write.csv(fin,paste0(output.dir, sampleName, "doubletAssignment.csv" ))
  rownames(fin) <- fin$cellID
  
  sample <- subset(sample, cells = barcodeCount$CellID)
  sample$doubletState <- fin$doubletState
  table(fin$doubletState)
  sample <- subset(sample, subset = doubletState == "singlet")
  ####################
  #Sample QC
  mito.pattern <- "^mt-"
    
  # calculate mito percentage
  mito.genes <- grep(pattern = mito.pattern, x = rownames(x = sample@assays$RNA), value = TRUE)
  sample[["percent.mito"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
  
  #plot QC metrics
  VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T)
    
  #subset based on mito count and feature cutoffs
  mitoCutoff = 20
  countsLower = 1000
  countsUpper = 100000
  featuresLower = 100
  featuresUpper = 7500
  sample <- subset(sample, cells = colnames(sample)[which(sample$percent.mito < mitoCutoff)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA > countsLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA < countsUpper)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA > featuresLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA < featuresUpper)])
  
  message("Filtered dimensions")
  print(dim(sample))
    
  VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T)
  FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mito")
  FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  sample <- SCTransform(sample, vars.to.regress = "percent.mito", verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use  = 10101)
  
  
  #Cell cycle scoring and regression based on that
  # perform cell cycle analysis
  
  sample <- NormalizeData(sample)
  sample <- CellCycleScoring(sample, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse, set.ident = TRUE, search = F)
    
  # add in cell cycle difference scoring
  sample$CC.Difference <- sample$S.Score - sample$G2M.Score
  
  # run PCA on cell cycle components
  sample <- RunPCA(sample, features = c(s.genes.mouse, g2m.genes.mouse))
  DimPlot(sample, reduction = "pca", group.by = 'Phase', label = T, label.size = 5)
  
  #Normalising to account for cell cycle and mitochondrial gene percentage
  regressVars = c("percent.mito", "S.Score", "G2M.Score")
  regressVars <- str_split(regressVars, ',')[[1]]
  
  
  sample <- SCTransform(sample, vars.to.regress = regressVars, verbose = T)
  
  
  #Dimensional reduction & clustering
    message("Running dimensional reduction algos")
    
    # always calculate 100 PCs
  sample <- RunPCA(object = sample, verbose = T, npcs = 100, seed.use = seed.use)
  
  # Elbow plot
  percent.var <- Stdev(sample)
  elbow.dim <- PCAtools::findElbowPoint(percent.var)
  print(elbow.dim)
    
  # find neighbors based on 1 to elbow dimensions
  sample <- FindNeighbors(object = sample, verbose = T, dims = 1:elbow.dim)
    
  # Run tSNE & UMAP on 1 to dims
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap"))
  
    
  # run a sequence of cluster resolutions for elbow dim
  sample <- FindClusters(object = sample, verbose = T, resolution = 0.8, random.seed = seed.use)
  
  
  ######################################
  #Cell Type Annotation for singlets
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  db = "~/Keerthana/CellTypeCountData/AdditionalData/SPLINTR_FinalData/SPLINTR-scTypeMarker.xlsx"
  tissue = "Leukemia"
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
  
  SaveSeuratRds(sample, paste0(outputRDS, sampleName, ".rds"))
  
  umap_coords <- Embeddings(sample, "umap")
  df <- data.frame(umap_coords)
  df$cluster <- sample$customclassif
  colnames(df) <- c("umap_1", "umap_2", "cluster")
  write.csv(df, paste0(outputDirectoryData, sampleName, ".csv"))
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
                    color = "black", bg.color = "white", bg.r = 0.1,
                    point.padding = unit(0.5, "lines")) +
    theme(legend.position = "none") +labs(title = NULL)
  
  plot
  
  ggsave(filename = paste0(outputDirectoryPlots, sampleName, ".png"), plot = plot, width = 6, height = 4)
  svglite(filename = paste0(outputDirectoryPlots, sampleName, ".svg"), width = 6, height = 4) #_withoutFiltering
  plot(plot)
  dev.off()
  
  cellType <- as.data.frame(table(df$cluster))
  write.csv(cellType, paste0(outputDirectoryData,sampleName,"_CellTypeCounts.csv") )
}
