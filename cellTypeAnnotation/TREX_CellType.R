lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(svglite)
#Path to directories
inputDirectory <- "~/Keerthana/CellTypeCountData/AdditionalData/TREX_FinalData/"
outputDirectory <- "~/Keerthana/CellTypeCountData/TREX/"
outputRDS <- "~/Keerthana/CellTypeCountData/finalRDS/TREX/"
outputDirectoryPlots <- "~/Keerthana/CellTypeCountData/plots/TREX/"
outputDirectoryData <- "~/Keerthana/CellTypeCountData/plotData/TREX/"

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
print(s.genes)
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
exclusiveMarkers <- c("Igf2", "Pf4", "Hexb", "Rsph1", "Pdgfra", "Bmp4", "Mog", "Clic6", "Rgs5", "Cldn5", "Reln", "Igfbpl1", "Slc32a1", "Slc17a7", "Aldoc")

sample = "brain1"
sampleList <- c("brain2", "brain3", "brain4")
for (i in sampleList){
  Processing_func(i)
}
#Preprocessing the TREX data
Processing_func <- function(sample){
  trex_input <- Read10X(data.dir = paste0(inputDirectory, "10X/", sample))
  # Initialize the Seurat object with the raw (non-normalized data).
  trex_brain <- CreateSeuratObject(counts = trex_input, min.cells = 3, min.features = 500)
  
  trex_brain <- subset(trex_brain, subset = nFeature_RNA > 500 & nFeature_RNA < 10000)
  
  ################################
  #Counting number of doublets and singlets
  singlets <- doubletFiltering(trex_brain)
  
  ######################################
  #Cell Type Annotation for singlets
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  db = "~/Keerthana/CellTypeCountData/AdditionalData/TREX/TREX_scType.xlsx"
  tissue = "Neural"
  gs_list = gene_sets_prepare(db, tissue)
  #singlet filtering
  trex_brain <- subset(trex_brain, cells = singlets$Cell.ID)
  
  trex_brain <- NormalizeData(trex_brain)
  trex_brain <- CellCycleScoring(trex_brain, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse, set.ident = TRUE)
  
  trex_brain$cc.difference = trex_brain$S.Score - trex_brain$G2M.Score
  trex_brain <- FindVariableFeatures(trex_brain, selection.method = "vst", nfeatures = 5000)
  trex_brain <- ScaleData(trex_brain, vars.to.regress = c("cc.difference"))
  trex_brain <- RunPCA(trex_brain)
  trex_brain <- FindNeighbors(trex_brain, dims = 1:20)
  trex_brain <- FindClusters(trex_brain, resolution = 0.8, n.start = 500, random.seed = 1024, algorithm = 1)
  trex_brain <- RunUMAP(trex_brain, dims = 1:20)
  es.max = sctype_score(scRNAseqData = trex_brain[["RNA"]]$scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = NULL) 
  
  cL_resutls = do.call("rbind", lapply(unique(trex_brain@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(trex_brain@meta.data[trex_brain@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(trex_brain@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  
  trex_brain@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    trex_brain@meta.data$customclassif[trex_brain@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  SaveSeuratRds(trex_brain, paste0(outputRDS, sample, ".rds"))
  
  umap_coords <- Embeddings(trex_brain, "umap")
  df <- data.frame(umap_coords)
  df$cluster <- trex_brain$customclassif
  df$cellID <- Cells(trex_brain)
  write.csv(df, paste0(outputDirectoryData, sample, ".csv"))
  write_csv(df, paste0("~/Keerthana/CellTypeCountData/plotData/UmapData/TREX/", sample))
  # colnames(df) <- c("cellID","umap_1", "umap_2", "cluster")
  # library(dplyr)
  # library(MASS)
  # library(ggrepel)
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
  # plot <- DimPlot(trex_brain, reduction = "umap", group.by = 'customclassif') + 
  #   geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
  #                   box.padding = unit(0.35, "lines"), fontface = "bold",
  #                   color = "black", bg.color = "white", bg.r = 0.1,
  #                   point.padding = unit(0.5, "lines")) +
  #   theme(legend.position = "none") +labs(title = NULL)
  # 
  # plot
  # 
  # ggsave(filename = paste0(outputDirectoryPlots, sample, ".png"), plot = plot, width = 6, height = 4)
  # svglite(filename = paste0(outputDirectoryPlots, sample, ".svg"), width = 6, height = 4) #_withoutFiltering
  # plot(plot)
  # dev.off()
  # 
  # cellType <- as.data.frame(table(df$cluster))
  # write.csv(cellType, paste0(outputDirectoryData,sample,"_CellTypeCounts.csv") )
  
}

#Doublet Filtering based on exclusive markers
doubletFiltering <- function(samplePaper){
  exclusiveMarkers <- c("Igf2", "Pf4", "Hexb", "Rsph1", "Pdgfra", "Bmp4", "Mog", "Clic6", "Rgs5", "Cldn5", "Reln", "Igfbpl1", "Slc32a1", "Slc17a7", "Aldoc")
  temp <- subset(samplePaper, features = exclusiveMarkers)
  
  mat = as.data.frame(temp[["RNA"]]$counts)
  result <- t(apply(mat, 1, calculate_threshold_and_binarize))
  # Binarized matrix
  binarized_matrix <- t(sapply(result, function(x) x[[1]]))
  
  # Stats matrix with row names
  stats_matrix <- t(sapply(result, function(x) unlist(x[[2]])))
  rownames(stats_matrix) <- rownames(mat)
  
  #Finding the number of genes expressed
  column_sums <- colSums(binarized_matrix)
  
  # Get cell names where 0, 1 and more than 1 exclusive marker is expressed
  doublets <- data.frame("Cell ID" = colnames(mat) [column_sums > 1])
  singlets_1 <- colnames(mat) [column_sums == 1]
  singlets_2 <- colnames(mat) [column_sums == 0]
  singlets <- data.frame("Cell ID" = c(singlets_1, singlets_2))
  return(singlets)
}

calculate_threshold_and_binarize <- function(row) {
  threshold = 0
  binarized_row <- ifelse(row > threshold, 1, 0)
  quartile_stats <- c(min = min(row), Q1 = quantile(row, probs = 0.25), median = median(row), mean = mean(row), Q3 = quantile(row, probs = 0.75), max = max(row), threshold = threshold)
  return(list(binarized_row, quartile_stats))
}
