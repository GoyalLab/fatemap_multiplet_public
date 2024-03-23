### generating doublets for Zhang, Melzer, et al., 2024 in order to make classifier
### Created by Madeline E Melzer on 20231211
### Last edited by Madeline E Melzer on 20240129
### code adapted from Ziyang Zhang, count_doublets_utils.py (chuckzzzz/doublet_benchmark github repo) and Goyal et al 2023, FM01_s1s2_integration_scTransForm_20210427_V2.R

library(tidyverse)
library(Seurat)
library(DropletUtils)
library(glue)
library(sctransform)
library(umap)
options(future.globals.maxSize = 4000 * 1024^2)

set.seed = 23

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
  assay_for_simulated_doublets <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")[, c(real_cells_pt1, real_cells_pt2)]
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


datasets_samples <- list(
  #"LARRY" = "d2_1", (not in yet)
  #"smartseq3_reads" = "brain1",
  #"smartseq3_umis" = "brain1",
  "SPLINTR" = "chemoDay2_1" #(no cell names (colnames) present in the input matrix)
)

dataset = "FM01"
sample = "sample1"

doublet_rates = list(0.05, 0.08, 0.15, 0.20, 0.25)

for (doublet_rate in doublet_rates) {
  
  data_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/10X/", sample, "/")
  singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/")
  output_dir <- paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/variable_doublet_rates_2/", doublet_rate)
  labeltable_dir <- paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/")
  
  datacreatedoublets = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = doublet_rate)
  data = PercentageFeatureSet(datacreatedoublets, pattern = "^MT-", col.name = "percent.mt")
  data_norm = NormalizeData(data)
  data_norm_scaled = ScaleData(data_norm, features = rownames(data_norm), vars.to.regress = "percent.mt")
  data_norm_scaled_var = FindVariableFeatures(data_norm_scaled)
  
  # save
  matrix <- data_norm_scaled_var@assays$RNA@counts
  features <- rownames(matrix) #note: swapped rows and columns
  barcodes <- colnames(matrix)
  labels = data_norm_scaled_var@meta.data$label
  
  Matrix::writeMM(matrix, file = paste0(output_dir, "matrix.mtx"))
  write.table(features, file = gzfile(paste0(output_dir, "features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(barcodes, file = gzfile(paste0(output_dir, "barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  metadata <- data.frame(barcode = barcodes, label = labels)
  write.csv(metadata, file = paste0(labeltable_dir, glue("labels_{doublet_rate}_2.csv")), row.names = FALSE)
}

dataset_samples <- list(
  #"FM01" = "sample1",
  #"FM02" = "1-1uMPLX",
  #"FM03"= "DMSO-FM3-1uMPLX",
  #"FM04" = "BC18_B1",
  #"FM05" = "FM05-250nMPLX-A1",
  #"FM06"= "FM06-WM989Naive-1",
  #"FM08" = "run1_sample3",
  #"non_cancer" = "10kbarsiblingA",
  #"s1nc_positiveControl",
  "Biorxiv" = "1_DMSO_A",
  #"TREX" = "brain1",
  #"SPLINTR" = "chemoDay2_1",
  "ClonMapper" = "FM1",
  "LARRY" = "LK_d2",
  #"cellTag" = "d2-RNA-5",
  #"watermelon" = "T47D-naive-1"
  #"smartseq3_reads" = "brain1",
  #"smartseq3_umis" = "brain1",
  #"SPLINTR" = "chemoDay2_1" #(no cell names (colnames) present in the input matrix)
)

doublet_rate = 0.1

for (dataset in names(dataset_samples)) {
  sample = dataset_samples[[dataset]]
  
  data_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/10X/", sample)
  singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/")
  output_dir <- paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/10X_doublets_2/")
  labeltable_dir <- paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/")
  
  datacreatedoublets = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = doublet_rate)
  data = PercentageFeatureSet(datacreatedoublets, pattern = "^MT-", col.name = "percent.mt")
  data_norm = NormalizeData(data)
  data_norm_scaled = ScaleData(data_norm, features = rownames(data_norm), vars.to.regress = "percent.mt")
  data_norm_scaled_var = FindVariableFeatures(data_norm_scaled)
  
  # save
  matrix <- data_norm_scaled_var@assays$RNA@counts
  features <- rownames(matrix) #note: swapped rows and columns
  barcodes <- colnames(matrix)
  labels = data_norm_scaled_var@meta.data$label
  
  Matrix::writeMM(matrix, file = paste0(output_dir, "matrix.mtx"))
  write.table(features, file = gzfile(paste0(output_dir, "features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(barcodes, file = gzfile(paste0(output_dir, "barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  metadata <- data.frame(barcode = barcodes, label = labels)
  write.csv(metadata, file = paste0(labeltable_dir, glue("labels_2.csv")), row.names = FALSE)
}









########## positive control: merging and scTransforming sample1 with non_cancer A in order to get positive control for training the classifier!

sample1_dir = "/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/FM01/10X/sample1"
non_cancer_dir = "/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/non_cancer/10X/10kbarsiblingA"

dataset = "FM01"
sample = "sample1"
setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/")
data_sample1 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.1)
singlets_sample1 <- subset(data_sample1, subset = label == "singlet")

dataset = "non_cancer"
sample = "10kbarsiblingA"
setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/")
data_non_cancer = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.01)
singlets_non_cancer <- subset(data_non_cancer, subset = label == "singlet")
singlets_non_cancer@meta.data$label <- "non_cancer"
selected_cells = sample(colnames(singlets_non_cancer), size = 0.1*length(colnames(singlets_sample1))) #subsampling to get 10% of the amount of cells as the singlets in FM01 sample1
singlets_non_cancer_subset <- subset(singlets_non_cancer, cells = selected_cells)


sample1_singlets = PercentageFeatureSet(singlets_sample1, pattern = "^MT-", col.name = "percent.mt")
non_cancer_singlets = PercentageFeatureSet(singlets_non_cancer_subset, pattern = "^MT-", col.name = "percent.mt")

s1nc_scTransform <- merge(sample1_singlets, y = non_cancer_singlets, add.cell.ids = c("S1", "NC"), project = "S1NC")

s1nc_scTransform.list <- SplitObject(s1nc_scTransform, split.by = "label")
s1nc_scTransform.list <- s1nc_scTransform.list[c("singlet", "non_cancer")]

for (i in 1:length(s1nc_scTransform.list)) {
  s1nc_scTransform.list[[i]] <- SCTransform(s1nc_scTransform.list[[i]])
}

s1nc_scTransform.features <- SelectIntegrationFeatures(object.list = s1nc_scTransform.list, nfeatures = 7000)
options(future.globals.maxSize = 1500 * 1024^2)
s1nc_scTransform.list <- PrepSCTIntegration(object.list = s1nc_scTransform.list, anchor.features = s1nc_scTransform.features, 
                                            verbose = FALSE)

s1nc_scTransform.anchors <- FindIntegrationAnchors(object.list = s1nc_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1nc_scTransform.features, verbose = FALSE)
s1nc_scTransform.integrated <- IntegrateData(anchorset = s1nc_scTransform.anchors, normalization.method = "SCT", verbose = FALSE)

s1nc_scTransform.integrated <- RunPCA(s1nc_scTransform.integrated, verbose = FALSE)
s1nc_scTransform.integrated <- RunUMAP(s1nc_scTransform.integrated, dims = 1:50)

DimPlot(s1nc_scTransform.integrated, reduction = "pca", group.by = "label", pt.size = 2) # UMAP plot colored by 'label'
DimPlot(s1nc_scTransform.integrated, reduction = "umap", group.by = "label", pt.size = 2)

#saveRDS(s1nc_scTransform.integrated, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1nc_positiveControl/s1nc_scTransform_50PCs_2.rds")

s1nc_scTransform.integrated = readRDS("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1nc_positiveControl/s1nc_scTransform_50PCs_2.rds")



### saving
setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1nc_positiveControl/10X_doublets_3/")
matrix <- s1nc_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1nc_scTransform.integrated@meta.data$label

Matrix::writeMM(matrix, file = "matrix.mtx")
#system("gzip matrix.mtx")
write.table(features, file = gzfile("features.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile("barcodes.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)


setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1nc_positiveControl/")
metadata <- data.frame(barcode = barcodes, label = labels)
write.csv(metadata, file = "labels_3.csv", row.names = FALSE)


#######################################

### sanity check for average transcriptomes 

data = PercentageFeatureSet(data_sample1, pattern = "^MT-", col.name = "percent.mt")
data_norm = NormalizeData(data)
data_norm_scaled = ScaleData(data_norm, features = rownames(data_norm), vars.to.regress = "percent.mt")
data_norm_scaled_var = FindVariableFeatures(data_norm_scaled)
#data = RunPCA(data_norm_scaled_var, verbose = FALSE)
#data = RunUMAP(data, features)

data.matrix = data@assays$RNA@counts

data.matrix.transposed = t(data.matrix)


data.umap.t = umap(data.matrix.transposed)

data.umap.df = data.umap.t$layout %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2"))

metadata <- data@meta.data
data.umap.df$label <- metadata[rownames(data.umap.df), "label"]
plot <- ggplot(data.umap.df, aes(x = UMAP1, y = UMAP2, color = label)) +
  geom_point(size = 3) +
  theme_classic(base_size = 18)

plot


#the wrong plot
data.umap.genes = umap(data.matrix)
data.umap.genes.df = data.umap.genes$layout %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2"))

data.umap.genes.df$label <- rownames(data[["RNA"]])
data.umap.genes.df$rowSums = rowSums(data.matrix)
filtered_data <- filter(data.umap.genes.df, rowSums > 6000 & rowSums < 20000)

plot <- ggplot(filtered_data, aes(x = UMAP1, y = UMAP2, color = rowSums)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_classic(base_size = 18)

plot








gene_of_interest <- "TAGLN"  # replace GENE_NAME with your gene name


data.umap.genes.df$color <- ifelse(data.umap.genes.df$label == gene_of_interest, "red", "black")

plot <- ggplot() +
  theme_classic(base_size = 18)

plot <- plot + geom_point(data = subset(data.umap.genes.df, label != gene_of_interest),
                          aes(x = UMAP1, y = UMAP2), color = "black")

# Add the points for label gene_of_interest on top
plot <- plot + geom_point(data = subset(data.umap.genes.df, label == gene_of_interest),
                          aes(x = UMAP1, y = UMAP2), color = "red")

plot

write.csv(data.umap.df, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/FM01/UMAP_cells.csv", row.names = FALSE)

### on DE genes

de = read.csv("/Users/mem3579/quest/ZhangMelzerEtAl/data/DE/FM01/control_3_vs._7__MAST.csv")


data.umap.genes.df = data.umap.genes.df %>% rename("label" = "gene")
data_genes_de = left_join(data.umap.genes.df, de, by = "gene")
left_join()

library(viridis)
plot <- ggplot(data_genes_de, aes(x = UMAP1, y = UMAP2, color = avg_log2FC)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_classic(base_size = 18)

plot


plot <- ggplot() +
  geom_point(data = subset(data_genes_de, is.na(avg_log2FC)), 
             aes(x = UMAP1, y = UMAP2), color = "grey", size = 3) +
  theme_classic(base_size = 18)

# Add points with existing avg_log2FC values on top
plot <- plot + geom_point(data = subset(data_genes_de, !is.na(avg_log2FC)), 
                          aes(x = UMAP1, y = UMAP2, color = avg_log2FC), size = 3) +
  scale_color_viridis_c()

plot




































############################### merging and scTransforming sample1 with sample2 to check fitting effect

sample1_dir = "/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/FM01/10X/sample1"
sample2_dir = "/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/FM01/10X/sample2"

dataset = "FM01"
sample = "sample1"
setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/")
sample1 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.1)
sample1@meta.data$sample <- "sample1"

sample = "sample2"
setwd(paste0("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/", dataset, "/"))
data_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/10X/", sample, "/")
singlet_list_dir <- paste0("/Users/mem3579/quest/ZhangMelzerEtAl/data/datasets/zz_all/", dataset, "/fatemapID/")
sample2 = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = 0.1)
sample2@meta.data$sample <- "sample2"

saveRDS(sample1, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/sample1/sample1.rds")
saveRDS(sample2, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/sample2/sample2.rds")

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

saveRDS(s1s2_scTransform.integrated, file = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/s1s2_scTransform_50PCs.rds")


# saving
setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/10X_doublets_2/")
matrix <- s1s2_scTransform.integrated@assays$integrated@data
features <- rownames(matrix)
barcodes <- colnames(matrix)
labels = s1s2_scTransform.integrated@meta.data$label
samples = s1s2_scTransform.integrated@meta.data$sample

Matrix::writeMM(matrix, file = "matrix.mtx")
#system("gzip matrix.mtx")
write.table(features, file = gzfile("features.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(barcodes, file = gzfile("barcodes.tsv.gz"), sep = "\t", row.names = FALSE, col.names = FALSE)


setwd("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/")
metadata <- data.frame(barcode = barcodes, label = labels, sample = samples)
write.csv(metadata, file = "labels_2.csv", row.names = FALSE)


# saving separate matrices for each sample 
metadata <- data.frame(barcode = colnames(s1s2_scTransform.integrated@assays$integrated@data),
                       label = s1s2_scTransform.integrated@meta.data$label,
                       sample = s1s2_scTransform.integrated@meta.data$sample)

samples_list <- split(metadata, metadata$sample)

for (sample_name in names(samples_list)) {
  # Filter the matrix and barcodes for the current sample
  current_sample_barcodes <- samples_list[[sample_name]]$barcode
  current_sample_matrix <- s1s2_scTransform.integrated@assays$integrated@data[, current_sample_barcodes]
  
  Matrix::writeMM(current_sample_matrix, file = paste0("matrix_", sample_name, ".mtx"))
  write.table(rownames(current_sample_matrix), file = gzfile(paste0("features_", sample_name, ".tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(current_sample_barcodes, file = gzfile(paste0("barcodes_", sample_name, ".tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Save the metadata for the current sample
  write.csv(samples_list[[sample_name]], file = paste0("labels_", sample_name, ".csv"), row.names = FALSE)
}

############ only getting the genes expressed in both samples
s1s2_scTransform.integrated = readRDS("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/s1s2/s1s2_scTransform_50PCs.rds")


counts_sample1 <- GetAssayData(sample1, slot = "counts")
counts_sample2 <- GetAssayData(sample2, slot = "counts")

total_counts_sample1 <- Matrix::rowSums(counts_sample1)
total_counts_sample2 <- Matrix::rowSums(counts_sample2)

nonzero_genes_sample1 <- names(total_counts_sample1[total_counts_sample1 > 0])
nonzero_genes_sample2 <- names(total_counts_sample2[total_counts_sample2 > 0])

overlapping_genes <- intersect(nonzero_genes_sample1, nonzero_genes_sample2)

s1s2_filtered <- subset(s1s2_scTransform.integrated, features = overlapping_genes)










##########################################################################################################################
# re-running high cross-sample cellID datasets
##########################################################################################################################

### for these, I need to change the function to make sure I am just getting the sample singlets. 

dataset_samples <- list(
  "Biorxiv" = "1_DMSO_B",
  #"ClonMapper" = "FM1",
  "LARRY" = "LK_d4_1"
  #"TREX" = "brain2"
)

doublet_rate = 0.1

for (dataset in names(dataset_samples)) {
  sample = dataset_samples[[dataset]]
  
  data_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/", dataset, "/10X/", sample, "/")
  singlet_list_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/datasets/zz_all/", dataset, "/fatemapID/")
  output_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/", dataset, "/10X_doublets_2/")
  labeltable_dir <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/", dataset, "/")
  
  datacreatedoublets = load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix = dataset, doublets_pct = doublet_rate)
  data = PercentageFeatureSet(datacreatedoublets, pattern = "^MT-", col.name = "percent.mt")
  data_norm = NormalizeData(data)
  data_norm_scaled = ScaleData(data_norm, features = rownames(data_norm), vars.to.regress = "percent.mt")
  data_norm_scaled_var = FindVariableFeatures(data_norm_scaled)
  
  # save
  matrix <- GetAssayData(data_norm_scaled_var, assay = "RNA", slot = "counts")
  features <- rownames(matrix) #note: swapped rows and columns
  barcodes <- colnames(matrix)
  labels = data_norm_scaled_var@meta.data$label
  
  Matrix::writeMM(matrix, file = paste0(output_dir, "matrix.mtx"))
  write.table(features, file = gzfile(paste0(output_dir, "features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(barcodes, file = gzfile(paste0(output_dir, "barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  metadata <- data.frame(barcode = barcodes, label = labels)
  write.csv(metadata, file = paste0(labeltable_dir, glue("labels_2.csv")), row.names = FALSE)
}


















