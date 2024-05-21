### checking crosss-experimental classifier results
### Created by Madeline E Melzer on 20240509
### Last edited by Madeline E Melzer on 20240515 (on Quest)
### code adapted (an is an extension of) from generatingDoublets.R, created by me.

library(tidyverse)
library(Seurat)
library(DropletUtils)
library(glue)
library(sctransform)
library(umap)
library(Matrix)
library(purrr)
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(23) #note, set 20240515

load_fatemap_into_seurat<-function(data_dir, cell_labels_prefix, doublets_pct){
  expression_matrix <- Read10X(data.dir = data_dir)
  
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  singlets_file<-glue("{data_dir}singlets_all.txt")
  
  # key file for simulating doublets
  singlets_pair_for_doublet_simulation_file<-list.files(data_dir, pattern="corrected_singlet_pairs.csv")
  
  # extract the singlets
  singlets_pair_for_doublet_simulation_file<-glue("{data_dir}/{singlets_pair_for_doublet_simulation_file}")
  # check if "-1" is a part of 10X barcodes
  if(grepl("-1", seurat_object|>colnames()|>head(1))){
    singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
      as.vector()|>
      paste("-1", sep = "")
  }else{
    singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
      as.vector()|>
      str_remove("-1")
  }
  
  # extract the singlet pairs
  singlet_pair_df<-read.csv(singlets_pair_for_doublet_simulation_file, header = F)|>
    #dplyr::mutate_each(funs = function(x) str_remove_all(x," "))
    dplyr::mutate(across(everything(), ~str_remove_all(.x, " ")))
  # remove combinations that did not pass QC in seurat
  all_seurat_cells<-seurat_object@meta.data|>row.names()
  
  # take average of expression to produce doublets
  real_cells_pt1<-singlet_pair_df$V1
  real_cells_pt2<-singlet_pair_df$V2
  #assay_for_simulated_doublets<-seurat_object@assays$RNA@counts[,c(real_cells_pt1, real_cells_pt2)]
  assay_for_simulated_doublets <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")[, c(real_cells_pt1, real_cells_pt2)]
  avg_assay<-(assay_for_simulated_doublets[, real_cells_pt1] + assay_for_simulated_doublets[, real_cells_pt2])/2
  
  # create artificial doublet ID
  singlets_in_seurat <- singlets[singlets %in% colnames(seurat_object)]
  doublet_id<-paste(singlet_pair_df$V1,singlet_pair_df$V2, sep = "--")
  colnames(avg_assay)<-doublet_id
  num_doublets<-round(length(singlets_in_seurat)*doublets_pct/(1+doublets_pct)) #derived from the case where 2 unique singlets make up 1 doublet, num_doublets = doublets_pct*(num_singlets-num_doublets)
  avg_assay<-avg_assay[,1:num_doublets]
  doublet_object<-CreateSeuratObject(counts = avg_assay)
  
  # remove singlets used to generate doublets
  singlets_used_for_doublets <- unique(unlist(strsplit(colnames(avg_assay), "--")))
  singlets_to_keep <- setdiff(colnames(seurat_object), singlets_used_for_doublets)
  singlets_to_keep_in_seurat <- intersect(singlets_in_seurat, singlets_to_keep)
  if (length(singlets_to_keep_in_seurat) > (num_doublets*(1-doublets_pct)/doublets_pct)) {
    singlets_trimmed = sample(singlets_to_keep_in_seurat, num_doublets*(1-doublets_pct)/doublets_pct)
  } else {
    num_doublets = round(length(singlets_in_seurat)*doublets_pct/(1-doublets_pct)) #calculate the number of doublets, no assumptions about removal being made
    avg_assay_new<-avg_assay[,1:num_doublets] 
    singlets_used_for_doublets <- unique(unlist(strsplit(colnames(avg_assay_new), "--")))
    singlets_to_keep <- setdiff(colnames(seurat_object), singlets_used_for_doublets) 
    singlets_to_keep_in_seurat <- intersect(singlets_in_seurat, singlets_to_keep) # remove singlets used to make doublets
    num_doublets_adjusted = round((length(singlets_to_keep_in_seurat)*doublets_pct) / (1-doublets_pct)) # adjust the number of doublets based on the number of singlets removed
    avg_assay_adjusted = avg_assay_new[,1:num_doublets_adjusted] # subset the avg_assay_new to be as long as the adjusted doublet count
    doublet_object<-CreateSeuratObject(counts = avg_assay_adjusted) #create doublet object
    singlets_trimmed = singlets_to_keep_in_seurat #no more singlets removed after the ones used to make doublets are removed
  }
  singlet_object <- seurat_object[, singlets_trimmed]
  
  # tag doublet labels back to the doublet object
  label<-rep("doublet",ncol(doublet_object))|>
    as.data.frame()
  rownames(label)<-colnames(doublet_object)
  doublet_object<-AddMetaData(doublet_object, label, col.name = "label")
  
  # merge the two objects
  label<-rep("singlet",ncol(singlet_object))|>
    as.data.frame()
  rownames(label)<-colnames(singlet_object)
  singlet_object<-AddMetaData(singlet_object, label, col.name = "label")
  
  singlet_countsMatrix <- GetAssayData(singlet_object, assay = "RNA", slot = "counts")
  doublet_countsMatrix <- GetAssayData(doublet_object, assay = "RNA", slot = "counts")
  combined_counts <- cbind(singlet_countsMatrix, doublet_countsMatrix)
  doublet.seurat <- CreateSeuratObject(counts = combined_counts, project = "Combined")
  
  singlet_metadata <- singlet_object@meta.data
  doublet_metadata <- doublet_object@meta.data
  combined_metadata <- rbind(singlet_metadata, doublet_metadata)
  doublet.seurat@meta.data <- combined_metadata
  
  #doublet.seurat<-merge(singlet_object, y = doublet_object, project = "combined_doublet_singlet")
  doublet_count <- sum(doublet.seurat@meta.data$label == "doublet")
  print(doublet_count)
  singlet_count <- sum(doublet.seurat@meta.data$label == "singlet")
  print(singlet_count)
  doublet.seurat
  
}

#outputDataDir = "/projects/p31666/melzer/ZhangMelzerEtAl/classifier/crossExperiment/"
outputDataDir = "/projects/b1042/GoyalLab/melzer/ZhangMelzerEtAl/classifier/crossExperiment/"

############################### merging and scTransforming TREX brain 1 and brain 2 to check fitting effect (i.e. cross-experiment) for different conditions ################ 20240509
#### this is run on Quest

TREXDir = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/TREX/10X/"
#sample1_name = "brain2"
#sample2_name = "brain3"


SPLINTRDir = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/SPLINTR/10X/"
#sample1_name = "chemoDay2_1"
#sample2_name = "chemoDay2_2"

ClonMapperDir = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/ClonMapper/10X/"
sample1_name = "FM1"
sample2_name = "FM7"


#dataset = "TREX"
#dataset = "SPLINTR"
dataset = "ClonMapper"

#setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
#data_dir <- paste0(TREXdir, "sample1", "/")
sample1_dir = glue(ClonMapperDir, sample1_name, "/")
#singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/") #dont need because singlet list is in the sample directory itself
sample1 = load_fatemap_into_seurat(data_dir = sample1_dir, cell_labels_prefix = dataset, doublets_pct = 0.1)
sample1@meta.data$dataset <- dataset
sample1@meta.data$sample <- "sample1"
#saveRDS(sample1, file = paste0(outputDataDir, dataset, "/", sample1_name, ".rds"))

#setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
#data_dir <- paste0(TREXDir, "sample2", "/")
sample2_dir = glue(ClonMapperDir, sample2_name, "/")
#singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/")
sample2 = load_fatemap_into_seurat(data_dir = sample2_dir, cell_labels_prefix = dataset, doublets_pct = 0.1)
sample2@meta.data$dataset <- dataset
sample2@meta.data$sample <- "sample2"
#saveRDS(sample2, file = paste0(outputDataDir, dataset, "/", sample2_name, ".rds"))


#scTransforming
sample1 = PercentageFeatureSet(sample1, pattern = "^MT-", col.name = "percent.mt")
sample2 = PercentageFeatureSet(sample2, pattern = "^MT-", col.name = "percent.mt")

s1s2_scTransform <- merge(sample1, y = sample2, add.cell.ids = c("S1", "S2"), project = "S1S2")

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "sample")
s1s2_scTransform.list <- s1s2_scTransform.list[c("sample1", "sample2")]

for (i in 1:length(s1s2_scTransform.list)) {
  s1s2_scTransform.list[[i]] <- SCTransform(s1s2_scTransform.list[[i]])
}

s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 7000)

s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, 
                                            verbose = FALSE)

s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1s2_scTransform.features, verbose = FALSE)
s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", verbose = FALSE)

s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

DimPlot(s1s2_scTransform.integrated, reduction = "pca", group.by = "label", pt.size = 2) # UMAP plot colored by 'label'
DimPlot(s1s2_scTransform.integrated, reduction = "umap", group.by = "label", pt.size = 2)

#saveRDS(s1s2_scTransform.integrated, file = paste0(outputDataDir, dataset, glue("/{dataset}_integrated.rds")))


### saving
#setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/10X_doublets_2/")
if (!dir.exists(paste0(outputDataDir, dataset, "/10X/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/"))
}


if (!dir.exists(paste0(outputDataDir, dataset, "/10X/integrated/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/integrated/"))
}

matrix <- s1s2_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1s2_scTransform.integrated@meta.data$label
samples = s1s2_scTransform.integrated@meta.data$sample

Matrix::writeMM(matrix, file = paste0(outputDataDir, dataset, "/10X/integrated/", "matrix.mtx"))
#system("gzip matrix.mtx")
write.table(features, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated/", "features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated/", "barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)


#setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/")
metadata <- data.frame(barcode = barcodes, label = labels, sample = samples)
write.csv(metadata, file = paste0(outputDataDir, dataset, "/10X/integrated/", "labels_2.csv"), row.names = FALSE)


# saving separate matrices for each sample
metadata <- data.frame(barcode = colnames(s1s2_scTransform.integrated@assays$integrated@data),
                       label = s1s2_scTransform.integrated@meta.data$label,
                       sample = s1s2_scTransform.integrated@meta.data$sample)

samples_list <- split(metadata, metadata$sample)

for (sample_name in names(samples_list)) {
  
  if (!dir.exists(paste0(outputDataDir, dataset, "/10X/", sample_name))) {
    dir.create(paste0(outputDataDir, dataset, "/10X/", sample_name))
  }
  
  # Filter the matrix and barcodes for the current sample
  current_sample_barcodes <- samples_list[[sample_name]]$barcode
  current_sample_matrix <- s1s2_scTransform.integrated@assays$integrated@data[, current_sample_barcodes]
  
  Matrix::writeMM(current_sample_matrix, file = paste0(outputDataDir, dataset, "/10X/", sample_name, "/matrix.mtx"))
  write.table(rownames(current_sample_matrix), file = gzfile(paste0(outputDataDir, dataset, "/10X/", sample_name, "/features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(current_sample_barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/", sample_name, "/barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Save the metadata for the current sample
  write.csv(samples_list[[sample_name]], paste0(outputDataDir, dataset, "/10X/", sample_name, "/labels_2.csv"), row.names = FALSE)
}


######## integrating s1-s4 (both vehicle and chemo samples from SPLINTR mice) together

SPLINTRDir = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/SPLINTR/10X/"
dataset = "SPLINTR"
sample1_name = "chemoVehicle_1"
sample2_name = "chemoVehicle_2"
sample3_name = "chemoDay2_1"
sample4_name = "chemoDay2_2"

sample1 = readRDS(file = paste0(outputDataDir, dataset, "/", sample1_name, ".rds"))
sample2 = readRDS(file = paste0(outputDataDir, dataset, "/", sample2_name, ".rds"))
sample3 = readRDS(file = paste0(outputDataDir, dataset, "/", sample3_name, ".rds"))
sample4 = readRDS(file = paste0(outputDataDir, dataset, "/", sample4_name, ".rds"))


sample1 <- PercentageFeatureSet(sample1, pattern = "^MT-", col.name = "percent.mt")
sample2 <- PercentageFeatureSet(sample2, pattern = "^MT-", col.name = "percent.mt")
sample3 <- PercentageFeatureSet(sample3, pattern = "^MT-", col.name = "percent.mt")
sample4 <- PercentageFeatureSet(sample4, pattern = "^MT-", col.name = "percent.mt")


s1s2_scTransform <- merge(sample1, y = list(sample2, sample3, sample4), add.cell.ids = c("S1", "S2", "S3", "S4"), project = "S1S2")

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "sample")
s1s2_scTransform.list <- s1s2_scTransform.list[c("sample1", "sample2", "sample3", "sample4")]

for (i in 1:length(s1s2_scTransform.list)) {
  s1s2_scTransform.list[[i]] <- SCTransform(s1s2_scTransform.list[[i]])
}

s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 7000)

s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, 
                                            verbose = FALSE)

s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1s2_scTransform.features, verbose = FALSE)
s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", verbose = FALSE)

s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

DimPlot(s1s2_scTransform.integrated, reduction = "pca", group.by = "label", pt.size = 2) # UMAP plot colored by 'label'
DimPlot(s1s2_scTransform.integrated, reduction = "umap", group.by = "label", pt.size = 2)

saveRDS(s1s2_scTransform.integrated, file = paste0(outputDataDir, dataset, glue("/{dataset}_integrated_s1s2s3s4.rds")))


### saving
#setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/10X_doublets_2/")
if (!dir.exists(paste0(outputDataDir, dataset, "/10X/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/"))
}


if (!dir.exists(paste0(outputDataDir, dataset, "/10X/integrated_s1s2s3s4/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/integrated_s1s2s3s4/"))
}

matrix <- s1s2_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1s2_scTransform.integrated@meta.data$label
samples = s1s2_scTransform.integrated@meta.data$sample

Matrix::writeMM(matrix, file = paste0(outputDataDir, dataset, "/10X/integrated_s1s2s3s4/", "matrix.mtx"))
#system("gzip matrix.mtx")
write.table(features, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated_s1s2s3s4/", "features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated_s1s2s3s4/", "barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)


#setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/")
metadata <- data.frame(barcode = barcodes, label = labels, sample = samples)
write.csv(metadata, file = paste0(outputDataDir, dataset, "/10X/integrated_s1s2s3s4/", "labels_2.csv"), row.names = FALSE)


# saving separate matrices for each sample
metadata <- data.frame(barcode = colnames(s1s2_scTransform.integrated@assays$integrated@data),
                       label = s1s2_scTransform.integrated@meta.data$label,
                       sample = s1s2_scTransform.integrated@meta.data$sample)

samples_list <- split(metadata, metadata$sample)

for (sample_name in names(samples_list)) {
  
  if (!dir.exists(paste0(outputDataDir, dataset, "/10X/", sample_name, "_s1s2s3s4"))) {
    dir.create(paste0(outputDataDir, dataset, "/10X/", sample_name, "_s1s2s3s4"))
  }
  
  # Filter the matrix and barcodes for the current sample
  current_sample_barcodes <- samples_list[[sample_name]]$barcode
  current_sample_matrix <- s1s2_scTransform.integrated@assays$integrated@data[, current_sample_barcodes]
  
  Matrix::writeMM(current_sample_matrix, file = paste0(outputDataDir, dataset, "/10X/", sample_name, "_s1s2s3s4", "/matrix.mtx"))
  write.table(rownames(current_sample_matrix), file = gzfile(paste0(outputDataDir, dataset, "/10X/", sample_name, "_s1s2s3s4","/features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(current_sample_barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/", sample_name, "_s1s2s3s4","/barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Save the metadata for the current sample
  write.csv(samples_list[[sample_name]], paste0(outputDataDir, dataset, "/10X/", sample_name, "_s1s2s3s4","/labels_2.csv"), row.names = FALSE)
}



####### integrating FateMap, TREX, SPLINTR, and ClonMapper samples together (can't! mouse and human!). 

dataset = "FateMap"
FateMap1 = readRDS(file = paste0(outputDataDir, dataset, "/", "sample1", ".rds"))
FateMap2 = readRDS(file = paste0(outputDataDir, dataset, "/", "sample2", ".rds"))

dataset = "TREX"
TREX1 = readRDS(file = paste0(outputDataDir, dataset, "/", "brain2", ".rds"))
TREX2 = readRDS(file = paste0(outputDataDir, dataset, "/", "brain3", ".rds"))

dataset = "SPLINTR"
SPLINTR1 = readRDS(file = paste0(outputDataDir, dataset, "/", "chemoVehicle_1", ".rds"))
SPLINTR2 = readRDS(file = paste0(outputDataDir, dataset, "/", "chemoVehicle_2", ".rds"))
SPLINTR3 = readRDS(file = paste0(outputDataDir, dataset, "/", "chemoDay2_1", ".rds"))
SPLINTR4 = readRDS(file = paste0(outputDataDir, dataset, "/", "chemoDay2_2", ".rds"))

dataset = "ClonMapper"
ClonMapper1 = readRDS(file = paste0(outputDataDir, dataset, "/", "FM1", ".rds"))
ClonMapper2 = readRDS(file = paste0(outputDataDir, dataset, "/", "FM7", ".rds"))

FateMap1 <- PercentageFeatureSet(FateMap1, pattern = "^MT-", col.name = "percent.mt")
FateMap2 <- PercentageFeatureSet(FateMap2, pattern = "^MT-", col.name = "percent.mt")
TREX1 <- PercentageFeatureSet(TREX1, pattern = "^MT-", col.name = "percent.mt")
TREX2 <- PercentageFeatureSet(TREX2, pattern = "^MT-", col.name = "percent.mt")
SPLINTR1 <- PercentageFeatureSet(SPLINTR1, pattern = "^MT-", col.name = "percent.mt")
SPLINTR2 <- PercentageFeatureSet(SPLINTR2, pattern = "^MT-", col.name = "percent.mt")
SPLINTR3 <- PercentageFeatureSet(SPLINTR3, pattern = "^MT-", col.name = "percent.mt")
SPLINTR4 <- PercentageFeatureSet(SPLINTR4, pattern = "^MT-", col.name = "percent.mt")
ClonMapper1 <- PercentageFeatureSet(ClonMapper1, pattern = "^MT-", col.name = "percent.mt")
ClonMapper2 <- PercentageFeatureSet(ClonMapper2, pattern = "^MT-", col.name = "percent.mt")

# Add metadata 'origin' to each dataset
FateMap1 <- AddMetaData(FateMap1, 'FateMap1', col.name = 'origin')
FateMap2 <- AddMetaData(FateMap2, 'FateMap2', col.name = 'origin')
TREX1 <- AddMetaData(TREX1, 'TREX1', col.name = 'origin')
TREX2 <- AddMetaData(TREX2, 'TREX2', col.name = 'origin')
SPLINTR1 <- AddMetaData(SPLINTR1, 'SPLINTR1', col.name = 'origin')
SPLINTR2 <- AddMetaData(SPLINTR2, 'SPLINTR2', col.name = 'origin')
SPLINTR3 <- AddMetaData(SPLINTR3, 'SPLINTR3', col.name = 'origin')
SPLINTR4 <- AddMetaData(SPLINTR4, 'SPLINTR4', col.name = 'origin')
ClonMapper1 <- AddMetaData(ClonMapper1, 'ClonMapper1', col.name = 'origin')
ClonMapper2 <- AddMetaData(ClonMapper2, 'ClonMapper2', col.name = 'origin')


# Add metadata 'pair' to each dataset
FateMap1 <- AddMetaData(FateMap1, 'A', col.name = 'pair')
FateMap2 <- AddMetaData(FateMap2, 'B', col.name = 'pair')
TREX1 <- AddMetaData(TREX1, 'A', col.name = 'pair')
TREX2 <- AddMetaData(TREX2, 'B', col.name = 'pair')
SPLINTR1 <- AddMetaData(SPLINTR1, 'A', col.name = 'pair')
SPLINTR2 <- AddMetaData(SPLINTR2, 'B', col.name = 'pair')
SPLINTR3 <- AddMetaData(SPLINTR3, 'A', col.name = 'pair')
SPLINTR4 <- AddMetaData(SPLINTR4, 'B', col.name = 'pair')
ClonMapper1 <- AddMetaData(ClonMapper1, 'A', col.name = 'pair')
ClonMapper2 <- AddMetaData(ClonMapper2, 'B', col.name = 'pair')

dataset = "all"

# Merge all mouse samples into a single Seurat object
s1s2_scTransform <- merge(TREX1, y = list(TREX2, SPLINTR1, SPLINTR2, SPLINTR3, SPLINTR4), add.cell.ids = c("TREX1", "TREX2", "SPLINTR1", "SPLINTR2", "SPLINTR3", "SPLINTR4"), project = "S1S2")

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "origin")
s1s2_scTransform.list <- s1s2_scTransform.list[c("TREX1", "TREX2", "SPLINTR1", "SPLINTR2", "SPLINTR3", "SPLINTR4")]

for (i in 1:length(s1s2_scTransform.list)) {
  s1s2_scTransform.list[[i]] <- SCTransform(s1s2_scTransform.list[[i]])
}

s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 7000)

s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, 
                                            verbose = TRUE)

s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1s2_scTransform.features, verbose = FALSE)
s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", verbose = FALSE)

s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

DimPlot(s1s2_scTransform.integrated, reduction = "pca", group.by = "label", pt.size = 2) # UMAP plot colored by 'label'
DimPlot(s1s2_scTransform.integrated, reduction = "umap", group.by = "label", pt.size = 2)

saveRDS(s1s2_scTransform.integrated, file = paste0(outputDataDir, dataset, glue("/{dataset}_integrated_mouse.rds")))

s1s2_scTransform.integrated = readRDS(file = paste0(outputDataDir, dataset, glue("/{dataset}_integrated_mouse.rds")))

### saving
#setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/10X_doublets_2/")
if (!dir.exists(paste0(outputDataDir, dataset, "/10X/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/"))
}


if (!dir.exists(paste0(outputDataDir, dataset, "/10X/integrated_mouse/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/integrated_mouse/"))
}

matrix <- s1s2_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1s2_scTransform.integrated@meta.data$label
samples = s1s2_scTransform.integrated@meta.data$sample
origins = s1s2_scTransform.integrated@meta.data$origin
pairs = s1s2_scTransform.integrated@meta.data$pair

Matrix::writeMM(matrix, file = paste0(outputDataDir, dataset, "/10X/integrated_mouse/", "matrix.mtx"))
#system("gzip matrix.mtx")
write.table(features, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated_mouse/", "features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated_mouse/", "barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)


metadata <- data.frame(barcode = barcodes, label = labels, sample = samples, origin = origins, pair = pairs)
write.csv(metadata, file = paste0(outputDataDir, dataset, "/10X/integrated_mouse/", "labels_2.csv"), row.names = FALSE)



# saving a matrix for pair A
metadata <- data.frame(barcode = colnames(s1s2_scTransform.integrated@assays$integrated@data),
                       label = s1s2_scTransform.integrated@meta.data$label,
                       sample = s1s2_scTransform.integrated@meta.data$sample,
                       origin = s1s2_scTransform.integrated@meta.data$origin,
                       pair = s1s2_scTransform.integrated@meta.data$pair)

pairs_list <- split(metadata, metadata$pair)
pairs_list <- purrr::keep(pairs_list, ~ any(.x$pair == "A"))

for (pair_name in names(pairs_list)) {
  
  if (!dir.exists(paste0(outputDataDir, dataset, "/10X/", pair_name, "_mouse"))) {
    dir.create(paste0(outputDataDir, dataset, "/10X/", pair_name, "_mouse"))
  }
  
  # Filter the matrix and barcodes for the current sample
  current_sample_barcodes <- pairs_list[[pair_name]]$barcode
  current_sample_matrix <- s1s2_scTransform.integrated@assays$integrated@data[, current_sample_barcodes]
  
  Matrix::writeMM(current_sample_matrix, file = paste0(outputDataDir, dataset, "/10X/", pair_name, "_mouse", "/matrix.mtx"))
  write.table(rownames(current_sample_matrix), file = gzfile(paste0(outputDataDir, dataset, "/10X/", pair_name, "_mouse","/features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(current_sample_barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/", pair_name, "_mouse","/barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Save the metadata for the current sample
  write.csv(pairs_list[[pair_name]], paste0(outputDataDir, dataset, "/10X/", pair_name, "_mouse","/labels_2.csv"), row.names = FALSE)
}



# saving separate matrices for the origins in pair A
metadata <- data.frame(barcode = colnames(s1s2_scTransform.integrated@assays$integrated@data),
                       label = s1s2_scTransform.integrated@meta.data$label,
                       sample = s1s2_scTransform.integrated@meta.data$sample,
                       origin = s1s2_scTransform.integrated@meta.data$origin,
                       pair = s1s2_scTransform.integrated@meta.data$pair)

origins_list <- split(metadata, metadata$origin)
origins_list <- purrr::keep(origins_list, ~ any(.x$pair == "A"))

for (origin_name in names(origins_list)) {
  
  if (!dir.exists(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_mouse"))) {
    dir.create(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_mouse"))
  }
  
  # Filter the matrix and barcodes for the current sample
  current_sample_barcodes <- origins_list[[origin_name]]$barcode
  current_sample_matrix <- s1s2_scTransform.integrated@assays$integrated@data[, current_sample_barcodes]
  
  Matrix::writeMM(current_sample_matrix, file = paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_mouse", "/matrix.mtx"))
  write.table(rownames(current_sample_matrix), file = gzfile(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_mouse","/features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(current_sample_barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_mouse","/barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Save the metadata for the current sample
  write.csv(origins_list[[origin_name]], paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_mouse","/labels_2.csv"), row.names = FALSE)
}









##### Merge all human samples into a single Seurat object
s1s2_scTransform <- merge(FateMap1, y = list(FateMap2, ClonMapper1, ClonMapper2), add.cell.ids = c("FateMap1", "FateMap2", "ClonMapper1", "ClonMapper2"), project = "S1S2")

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "origin")
s1s2_scTransform.list <- s1s2_scTransform.list[c("FateMap1", "FateMap2", "ClonMapper1", "ClonMapper2")]

for (i in 1:length(s1s2_scTransform.list)) {
  s1s2_scTransform.list[[i]] <- SCTransform(s1s2_scTransform.list[[i]])
}

s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 7000)

s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, 
                                            verbose = TRUE)

s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1s2_scTransform.features, verbose = FALSE)
s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", verbose = FALSE)

s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

DimPlot(s1s2_scTransform.integrated, reduction = "pca", group.by = "label", pt.size = 2) # UMAP plot colored by 'label'
DimPlot(s1s2_scTransform.integrated, reduction = "umap", group.by = "label", pt.size = 2)

saveRDS(s1s2_scTransform.integrated, file = paste0(outputDataDir, dataset, glue("/{dataset}_integrated_human.rds")))

s1s2_scTransform.integrated = readRDS(file = paste0(outputDataDir, dataset, glue("/{dataset}_integrated_human.rds")))

### saving

if (!dir.exists(paste0(outputDataDir, dataset, "/10X/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/"))
}


if (!dir.exists(paste0(outputDataDir, dataset, "/10X/integrated_human/"))) {
  dir.create(paste0(outputDataDir, dataset, "/10X/integrated_human/"))
}

matrix <- s1s2_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1s2_scTransform.integrated@meta.data$label
samples = s1s2_scTransform.integrated@meta.data$sample

Matrix::writeMM(matrix, file = paste0(outputDataDir, dataset, "/10X/integrated_human/", "matrix.mtx"))
#system("gzip matrix.mtx")
write.table(features, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated_human/", "features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/integrated_human/", "barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)


metadata <- data.frame(barcode = barcodes, label = labels, sample = samples)
write.csv(metadata, file = paste0(outputDataDir, dataset, "/10X/integrated_human/", "labels_2.csv"), row.names = FALSE)




# saving a matrix for pair A
metadata <- data.frame(barcode = colnames(s1s2_scTransform.integrated@assays$integrated@data),
                       label = s1s2_scTransform.integrated@meta.data$label,
                       sample = s1s2_scTransform.integrated@meta.data$sample,
                       origin = s1s2_scTransform.integrated@meta.data$origin,
                       pair = s1s2_scTransform.integrated@meta.data$pair)

pairs_list <- split(metadata, metadata$pair)

for (pair_name in names(pairs_list)) {
  
  if (!dir.exists(paste0(outputDataDir, dataset, "/10X/", pair_name, "_human"))) {
    dir.create(paste0(outputDataDir, dataset, "/10X/", pair_name, "_human"))
  }
  
  # Filter the matrix and barcodes for the current sample
  current_sample_barcodes <- pairs_list[[pair_name]]$barcode
  current_sample_matrix <- s1s2_scTransform.integrated@assays$integrated@data[, current_sample_barcodes]
  
  Matrix::writeMM(current_sample_matrix, file = paste0(outputDataDir, dataset, "/10X/", pair_name, "_human", "/matrix.mtx"))
  write.table(rownames(current_sample_matrix), file = gzfile(paste0(outputDataDir, dataset, "/10X/", pair_name, "_human","/features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(current_sample_barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/", pair_name, "_human","/barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Save the metadata for the current sample
  write.csv(pairs_list[[pair_name]], paste0(outputDataDir, dataset, "/10X/", pair_name, "_human","/labels_2.csv"), row.names = FALSE)
}






# saving separate matrices for the origins in pair A
metadata <- data.frame(barcode = colnames(s1s2_scTransform.integrated@assays$integrated@data),
                       label = s1s2_scTransform.integrated@meta.data$label,
                       sample = s1s2_scTransform.integrated@meta.data$sample,
                       origin = s1s2_scTransform.integrated@meta.data$origin,
                       pair = s1s2_scTransform.integrated@meta.data$pair)

origins_list <- split(metadata, metadata$origin)
origins_list <- purrr::keep(origins_list, ~ any(.x$pair == "A"))

for (origin_name in names(origins_list)) {
  
  if (!dir.exists(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_human"))) {
    dir.create(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_human"))
  }
  
  # Filter the matrix and barcodes for the current sample
  current_sample_barcodes <- origins_list[[origin_name]]$barcode
  current_sample_matrix <- s1s2_scTransform.integrated@assays$integrated@data[, current_sample_barcodes]
  
  Matrix::writeMM(current_sample_matrix, file = paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_human", "/matrix.mtx"))
  write.table(rownames(current_sample_matrix), file = gzfile(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_human","/features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(current_sample_barcodes, file = gzfile(paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_human","/barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Save the metadata for the current sample
  write.csv(origins_list[[origin_name]], paste0(outputDataDir, dataset, "/10X/", origin_name, "_A_human","/labels_2.csv"), row.names = FALSE)
}




