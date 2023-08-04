set.seed(2022)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(PRROC)
library(pbapply)
library(scDblFinder)
library(data.table)
library(glue)
library(stringr)


cell_labels_prefix<-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM03/fatemapID/FM03"
sample_dir<-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM03/10X"
data_dir<-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM03"
samples <- c("DOT1Li-FM3-1uMPLX", "DMSO-FM3-1uMPLX")
all_sample_files <- glue("{sample_dir}/{samples}")

load_fatemap_into_seurat <- function(data_dir, cell_labels_prefix, doublets_pct = 0.08) {
  expression_matrix <- Read10X(data.dir = data_dir)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  singlets_file <- glue("{cell_labels_prefix}_singlets.txt")

  # key file for simulaitng doublets
  singlets_pair_for_doublet_simulation_file <- list.files(data_dir, pattern = "corrected_singlet_pairs.csv")

  # extract the singlets
  singlets_pair_for_doublet_simulation_file <- glue("{data_dir}/{singlets_pair_for_doublet_simulation_file}")
  singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
    as.vector()|>
    paste("-1", sep = "")

  # extract the singlet paris
  singlet_pair_df <- read.csv(singlets_pair_for_doublet_simulation_file, header = F)|>
    dplyr::mutate_each(funs = function(x) str_remove_all(x, " "))
  # remove combinations that did not pass QC in seurat
  all_seurat_cells <- seurat_object@meta.data|>row.names()

  # take average of expression to produce doublets
  real_cells_pt1 <- singlet_pair_df$V1
  real_cells_pt2 <- singlet_pair_df$V2
  assay_for_simulated_doublets <- seurat_object@assays$RNA@counts[, c(real_cells_pt1, real_cells_pt2)]
  avg_assay <- (assay_for_simulated_doublets[, real_cells_pt1] + assay_for_simulated_doublets[, real_cells_pt2]) / 2


  # create artificial doublet ID
  singlets_in_seurat <- singlets[singlets %in% colnames(seurat_object)]
  doublet_id <- paste(singlet_pair_df$V1, singlet_pair_df$V2, sep = "--")
  colnames(avg_assay) <- doublet_id
  num_doublets <- round(length(singlets_in_seurat) * doublets_pct / (1 - doublets_pct))
  avg_assay <- avg_assay[, 1:num_doublets]
  doublet_object <- CreateSeuratObject(counts = avg_assay)

  # tag doublet labels back to the doublet object
  label <- rep("doublet", ncol(doublet_object))|>
    as.data.frame()
  rownames(label) <- colnames(doublet_object)
  doublet_object <- AddMetaData(doublet_object, label, col.name = "label")

  # merge the two objects
  singlet_object <- seurat_object[, singlets_in_seurat]
  label <- rep("singlet", ncol(singlet_object))|>
    as.data.frame()
  rownames(label) <- colnames(singlet_object)
  singlet_object <- AddMetaData(singlet_object, label, col.name = "label")
  doublet.seurat <- merge(singlet_object, y = doublet_object, project = "combined_doublet_singlet")
  doublet.seurat
}

# FM03_s_list <- sapply(all_sample_files, function(x){
#   cur_s <- load_fatemap_into_seurat(x, cell_labels_prefix = cell_labels_prefix)
#   cur_s@meta.data[["sample"]] <- strsplit(x, "/")[[1]][10]
#   cur_s <- NormalizeData(cur_s)
#   cur_s <- FindVariableFeatures(cur_s, selection.method = "vst", nfeatures = 2000)
#   cur_s
# })
# features <- SelectIntegrationFeatures(object.list = FM03_s_list)
# anchors <- FindIntegrationAnchors(object.list = FM03_s_list, anchor.features = features)
# FM03_s <- IntegrateData(anchorset = anchors)
# DefaultAssay(FM03_s) <- "integrated"
# FM03_s <- ScaleData(FM03_s, verbose = FALSE)
# FM03_s <- RunPCA(FM03_s, npcs = 30, verbose = FALSE)
# FM03_s <- RunUMAP(FM03_s, reduction = "pca", dims = 1:30)
# FM03_s <- FindNeighbors(FM03_s, reduction = "pca", dims = 1:30)
# FM03_s <- FindClusters(FM03_s, resolution = 0.5)
# DefaultAssay(FM03_s)<-"RNA"
# saveRDS(FM03_s, "/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp/FM03_integrated.rds")


cell_labels_prefix<-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM04/fatemapID/FM04"
sample_dir<-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM04/10X"
data_dir<-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM04"
samples <- c("BC18_B2", "BC18_B1")
all_sample_files <- glue("{sample_dir}/{samples}")

FM04_s_list <- sapply(all_sample_files, function(x){
  cur_s <- load_fatemap_into_seurat(x, cell_labels_prefix = cell_labels_prefix)
  cur_s@meta.data[["sample"]] <- strsplit(x, "/")[[1]][10]
  cur_s <- NormalizeData(cur_s)
  cur_s <- FindVariableFeatures(cur_s, selection.method = "vst", nfeatures = 2000)
  cur_s
})
features <- SelectIntegrationFeatures(object.list = FM04_s_list)
anchors <- FindIntegrationAnchors(object.list = FM04_s_list, anchor.features = features)
FM04_s <- IntegrateData(anchorset = anchors)
DefaultAssay(FM04_s) <- "integrated"
FM04_s <- ScaleData(FM04_s, verbose = FALSE)
FM04_s <- RunPCA(FM04_s, npcs = 30, verbose = FALSE)
FM04_s <- RunUMAP(FM04_s, reduction = "pca", dims = 1:30)
FM04_s <- FindNeighbors(FM04_s, reduction = "pca", dims = 1:30)
FM04_s <- FindClusters(FM04_s, resolution = c(0.5,0.6,0.7,0.8))
DefaultAssay(FM04_s)<-"RNA"
saveRDS(FM04_s, "/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp/FM04_integrated.rds")
