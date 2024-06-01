#Neighbor Analysis for Barcode Multiplet Project
#Last updated 20240531 by Madeline E Melzer

### load libraries

library(tidyverse)
library(Seurat)
library(ggplot2)
library(reshape2)
library(svglite)
library(installr)
library(patchwork)
library(RANN)
library(DescTools)
library(tools)

### specify directories

setwd("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/My Drive/neighborAnalysis/")
dataDirectory = "/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/My Drive/neighborAnalysis/data/" #where your data files are
data2Directory = "/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/My Drive/neighborAnalysis/data/FM01/barcode/" #where your barcoded singlets and multiplet lists are
tableDirectory = "/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/My Drive/neighborAnalysis/data/tables/" #table output location
plotDirectory = "/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/My Drive/neighborAnalysis/plots/FM01/sample3/" #where your plots will go

set.seed(23)

# all elements that need to be renamed before saving are marked with "#### RENAME" 

#### Load the dataset

fm01.data <- Read10X(data.dir = "/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/My Drive/neighborAnalysis/data/FM01/10X/sample3/") #### RENAME

### Initialize the Seurat object with the raw (non-normalized data).
fm01 <- CreateSeuratObject(counts = fm01.data, project = "neighbor", min.cells = 3, min.features = 200)

### Quality Control
fm01[["percent.mt"]] <- PercentageFeatureSet(fm01, pattern = "^MT-")
VlnPlot(fm01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(fm01, feature1 = "nCount_RNA", feature2 = "percent.mt")

fm01 <- subset(fm01, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20) # thresholding values will be different for every dataset
fm01 <- NormalizeData(fm01)
fm01 <- FindVariableFeatures(fm01, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(fm01)
fm01 <- ScaleData(fm01, features = all.genes)
fm01 <- RunPCA(fm01, features = VariableFeatures(object = fm01))

### could examine and visualize PCA results a few different ways
    #print(fm01[["pca"]], dims = 1:5, nfeatures = 5)
    #VizDimLoadings(fm01, dims = 1:2, reduction = "pca")
    #DimPlot(fm01, reduction = "pca")
    #ElbowPlot(fm01)

### Get PCA coordinates
pcaCoordinates = fm01@reductions[['pca']]@cell.embeddings

cells_PCA = rownames(pcaCoordinates) #CellIds with Sample number as prefix
cells_PCA_cellID = sub("-1", "", cells_PCA) # Removing the -1 from cell ID
pcaCoordinates = as_tibble(pcaCoordinates) #UMAP coordinates in tibble
pcaCoordinates = pcaCoordinates %>% mutate(cellID = cells_PCA_cellID)
write.table(pcaCoordinates, file=paste0(dataDirectory,'pcaCoordinates.tsv'), col.names = TRUE, sep='\t')

### select only barcoded singlets

singletList = as_tibble(read.table(file = paste0(data2Directory, "FM01_singlets.txt"), header = FALSE, stringsAsFactors=F)) %>% rename(cellID = V1) #### RENAME
singletPCACoordinates = inner_join(pcaCoordinates, singletList, by = "cellID") #these cells are barcoded

### check out the PCA
ggplot(singletPCACoordinates, aes(x = PC_1, y = PC_2)) +
  geom_point() +
  labs(x = "PC_1", y = "PC_2") +
  ggtitle("Singlet Principal Components 1 vs 2")

### Generating average distance of 5NN for all barcoded singlets
singletPCA = as.data.frame(singletPCACoordinates[1:50])

averageDistances_singlets = c()
cvs_singlets = c()

for (i in 1:nrow(singletPCA)) {
  
  sampledCell = as.data.frame(singletPCA[i, , drop = FALSE])
  
  neighbors = as_tibble(nn2(singletPCA, query = sampledCell, k= 6))
  
  averageDistance = mean(neighbors$nn.dists[, 2:6])
  cv = CoefVar(as.matrix(neighbors$nn.dists[, 2:6]))
  
  averageDistances_singlets[i] = averageDistance
  cvs_singlets[i] = cv

  }

avgDistances_singlets = as_tibble(averageDistances_singlets)
cvs_singlets = as_tibble(cvs_singlets)

### 3 subsamples of 5NN distances for all cells (singlets and multiplets) where n = number of singlets 

allPCA <- as.data.frame(pcaCoordinates[1:50])
singletPCA <- as.data.frame(singletPCACoordinates[1:50])

averageDistances_allSims <- matrix(nrow = nrow(singletPCA), ncol = 3) #making a place to put the outputs for each "simulation"
cvs_allSims = matrix(nrow = nrow(singletPCA), ncol = 3)

for (i in 1:3) {
  randomIndices = sample(nrow(allPCA), nrow(singletPCA))
  randomCells = allPCA[randomIndices, ]
  
  averageDistances = c()
  cvs = c()
  
  for (j in 1:nrow(averageDistances_allSims)) {
    sampledCell = as.data.frame(randomCells[j, , drop = FALSE])
    
    neighbors = as_tibble(nn2(randomCells, query = sampledCell, k = 6))
    
    averageDistance = mean(neighbors$nn.dists[, 2:6])
    cv = CoefVar(as.matrix(neighbors$nn.dists[, 2:6]))
    
    averageDistances[j] = averageDistance
    cvs[j] = cv
    }
  
  averageDistances_allSims[ , i] = averageDistances
  cvs_allSims[ ,i] = cvs
  }

averageDistances_allSims = as.data.frame(averageDistances_allSims)
cvs_allSims = as.data.frame(cvs_allSims)


### Formatting for plotting

averageDistances_allSims_stack = sample_n(stack(averageDistances_allSims), nrow(avgDistances_singlets)) %>% select(-ind)
averageDistanceNorm <- mean(averageDistances_allSims_stack$values) #average of all cells subsampling (sims) to use for normalization

averageDistances_allSims_stack = averageDistances_allSims_stack %>% mutate(values = values - averageDistanceNorm)
avgDistances_singlets = avgDistances_singlets %>%  mutate(value = value - averageDistanceNorm)
averageDistances_partial = cbind(avgDistances_singlets, averageDistances_allSims_stack)

cvs_allSims_stack = sample_n(stack(cvs_allSims), nrow(cvs_singlets)) %>% select(-ind)
cvs_partial = cbind(cvs_singlets, cvs_allSims_stack)

typeA = c("singlets", "sim")

averageDistances_partialTrans = t(averageDistances_partial)
averageDistances_partialTitled = cbind(typeA, averageDistances_partialTrans)
averageDistances_partialTitled <- as.data.frame(averageDistances_partialTitled)
averageDistances_partialTitled$rowname <- rownames(averageDistances_partialTitled)

cvs_partialTrans = t(cvs_partial)
cvs_partialTitled = cbind(typeA, cvs_partialTrans)
cvs_partialTitled <- as.data.frame(cvs_partialTitled)
cvs_partialTitled$rowname <- rownames(cvs_partialTitled)

averageDistances_partial <- melt(averageDistances_partialTitled, id.vars = c("rowname", "typeA"), variable.name = "Column", value.name = "Value")
averageDistances_partial$rowname <- NULL
averageDistances_partial$Column <- NULL

cvs_partial <- melt(cvs_partialTitled, id.vars = c("rowname", "typeA"), variable.name = "Column", value.name = "Value")
cvs_partial$rowname <- NULL
cvs_partial$Column <- NULL

distances.summary = averageDistances_partial %>% group_by(typeA) %>% summarise(sd = sd(as.numeric(Value), na.rm = TRUE), Value = mean(as.numeric(Value), na.rm = TRUE))
cvs.summary = cvs_partial %>% group_by(typeA) %>% summarise(sd = sd(as.numeric(Value), na.rm = TRUE), Value = mean(as.numeric(Value), na.rm = TRUE))


plot = ggplot(averageDistances_partial, aes(x = typeA, y = as.numeric(Value))) +
  geom_point(aes(color = typeA), shape = 16, alpha = 0.7, size = 1, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 0.5)) +
  geom_violin(color = "black", alpha = 0, size = 0.8) +
  theme_classic((base_size = 18)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Subsampling Type", y = "Average Distance to 5NN") +
  geom_pointrange(aes(ymin = Value - sd, ymax = Value + sd, color = "black"), size = 1, color = "black", position = position_dodge(width = 0.5), data = distances.summary) + #keeping data.summary for non-subsampled points!
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'NJ01_filtered_feature_bc_matrix_6_distance.svg'), width = 4, height = 8) #### RENAME

plot = ggplot(cvs_partial, aes(x = typeA, y = as.numeric(Value))) +
  geom_point(aes(color = typeA), shape = 16, alpha = 0.7, size = 1, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 0.5)) +
  geom_violin(color = "black", alpha = 0, size = 0.8) +
  theme_classic((base_size = 18)) +
  labs(x = "Subsampling Type", y = "cv of Average Distance to 5NN") +
  geom_pointrange(aes(ymin = Value - sd, ymax = Value + sd, color = "black"), size = 1, color = "black", position = position_dodge(width = 0.5), data = cvs.summary) + #keeping data.summary for non-subsampled points!
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'NJ01_filtered_feature_bc_matrix_6_cv.svg'), width = 4, height = 8) #### RENAME

#Saving the table!

final = mutate(averageDistances_partial, sampleID = "NJ01_filtered_feature_bc_matrix_6") #### RENAME
write.table(final, file = paste0("/Users/mem3579/My Drive/neighborAnalysis/data/tables/", "NJ01_filtered_feature_bc_matrix_6.csv"), sep = ",", row.names = FALSE, col.names = TRUE) #### RENAME







###################### PLOTTING ALL DATASETS AND SAMPLES TOGETHER ###################################################################################################################################################################

plotDirectory = "/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/My Drive/neighborAnalysis/plots/" #### RENAME

csvFiles <- list.files(tableDirectory, pattern = ".csv", full.names = TRUE)

# Create an empty list to store the tables
tables_list <- list()

# Read each CSV file and store it in the list
for (file in csvFiles) {
  table <- read_csv(file)
  tables_list[[file]] <- table
}

# Combine all tables into a single dataframe
combined_df <- bind_rows(tables_list)

#they are already individually normalized to sim so you only have to plot singlets

neighbors.summary <- combined_df %>%
  mutate(sampleID = ifelse(str_starts(sampleID, "NJ01"), "NJ01_allsamples", sampleID)) %>%
  mutate(sampleID = ifelse(str_starts(sampleID, "FM08"), "FM08_allsamples", sampleID)) %>%
  group_by(sampleID, typeA) %>%
  summarise(
    sd = sd(as.numeric(Value), na.rm = TRUE),
    mean = mean(as.numeric(Value), na.rm = TRUE),
    median = median(as.numeric(Value), na.rm = TRUE),
    count = n()
  ) %>%
  filter(typeA == "singlets")
  

plot = ggplot(neighbors.summary, aes(x = sampleID, y = mean)) +
  geom_point(aes(color = sampleID), shape = 16, alpha = 0.7, size = 2) +
  #geom_violin(color = "black", alpha = 0, size = 0.8) +
  theme_classic(base_size = 18) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd, color = sampleID), shape = 16, size = 0.5) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(plot, file = paste0(plotDirectory, 'combinedMean.svg'), width = 8, height = 8)


plot = ggplot(neighbors.summary, aes(x = sampleID, y = median)) +
  geom_point(aes(color = sampleID), shape = 16, alpha = 0.7, size = 2) +
  theme_classic(base_size = 18) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  geom_pointrange(aes(ymin = median - sd, ymax = median + sd, color = sampleID), size = 0.5) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(plot, file = paste0(plotDirectory, 'combinedMedian.svg'), width = 8, height = 8)


###################### CONTROLS (only done for FM01_sample3) ###################################################################################################################################################################

pcaCoordinates_ordered = pcaCoordinates[order(pcaCoordinates$PC_1), ]

### looking at the PCA
ggplot(pcaCoordinates_ordered, aes(x = PC_1, y = PC_2)) +
  geom_point() +
  labs(x = "PC_1", y = "PC_2") +
  ggtitle("Singlet Principal Components 1 vs 2")


### taking 3 regions of the PCA plot and using them to create an unevenly distributed dataset
biasedSubsetA <- pcaCoordinates_ordered[pcaCoordinates_ordered$PC_1 >= -12.5 & pcaCoordinates_ordered$PC_1 <= -7.5 & pcaCoordinates_ordered$PC_2 >= -4 & pcaCoordinates_ordered$PC_2 <= -2, ]
biasedSubsetB <- pcaCoordinates_ordered[pcaCoordinates_ordered$PC_1 >= 2.5 & pcaCoordinates_ordered$PC_1 <= 7.5 & pcaCoordinates_ordered$PC_2 >= 8 & pcaCoordinates_ordered$PC_2 <= 10,  ]
biasedSubsetC <- pcaCoordinates_ordered[pcaCoordinates_ordered$PC_1 >= 12.5 & pcaCoordinates_ordered$PC_1 <= 17.5 & pcaCoordinates_ordered$PC_2 >= -11 & pcaCoordinates_ordered$PC_2 <= -9,  ]

pcaCoordinates_bias = rbind(biasedSubsetA, biasedSubsetB, biasedSubsetC)
pcaCoordinates_bias <- as.data.frame(pcaCoordinates_bias[1:50])

singletPCA <- as.data.frame(singletPCACoordinates[1:50]) #same as in previous section

### Determining neighbor distances for the biased subset

averageDistances_bias = c()
cvs_bias = c()

for (i in 1:nrow(pcaCoordinates_bias)) {
  
  sampledCell = as.data.frame(pcaCoordinates_bias[i, , drop = FALSE])
  neighbors = as_tibble(nn2(pcaCoordinates_bias, query = sampledCell, k= 6))
  
  averageDistance = mean(neighbors$nn.dists[, 2:6])
  cv = CoefVar(as.matrix(neighbors$nn.dists[, 2:6]))
  
  averageDistances_bias[i] = averageDistance
  cvs_bias[i] = cv
  
  }

avgDistances_bias = as_tibble(averageDistances_bias)
cvs_bias = as_tibble(cvs_bias)

avgDistances_bias$typeA = "bias"
cvs_bias$typeA = "bias"

avgDistances_singlets$typeA = "singlets"
cvs_singlets$typeA = "singlets"

averageDistances_total_bias = rbind(avgDistances_singlets, avgDistances_bias)
cvs_total_bias = rbind(cvs_singlets, cvs_bias)

### subsampling n = biased cells for the singlets because the biased subset n < singlet n

singletPCA_random = as.data.frame(sample_n(singletPCACoordinates[1:50], nrow(pcaCoordinates_bias)))

averageDistances_singletsSub = c()
cvs_singletsSub = c()

for (i in 1:nrow(singletPCA_random)) {
  
  sampledCell = as.data.frame(singletPCA_random[i, , drop = FALSE])
  neighbors = as_tibble(nn2(singletPCA_random, query = sampledCell, k= 6))
  
  averageDistance = mean(neighbors$nn.dists[, 2:6])
  cv = CoefVar(as.matrix(neighbors$nn.dists[, 2:6]))
  
  averageDistances_singletsSub[i] = averageDistance
  cvs_singletsSub[i] = cv
  
  }

avgDistances_singletsSub = as_tibble(averageDistances_singletsSub)
cvs_singletsSub = as_tibble(averageDistances_singletsSub)

avgDistances_singletsSub$typeA = "singlets"
cvs_singletsSub$typeA = "singlets"

### subsampling n = biased cells for all the cells (including singlets and multiplets) because the biased subset n < all cell n

allPCA <- as.data.frame(pcaCoordinates[1:50])

simPCA_random = as.data.frame(sample_n(pcaCoordinates[1:50], nrow(pcaCoordinates_bias)))

averageDistances_subSims = matrix(nrow = nrow(simPCA_random), ncol = 1) #making a place to put the outputs for each "simulation"
cvs_subSims = matrix(nrow = nrow(simPCA_random), ncol = 10)

for (i in 1) {
  randomIndices = sample(nrow(allPCA), nrow(simPCA_random))
  randomCells = allPCA[randomIndices, ]
  
  averageDistances = c()
  cvs = c()
  
  for (j in 1:nrow(averageDistances_subSims)) {
    sampledCell = as.data.frame(randomCells[j, , drop = FALSE])
    
    neighbors = as_tibble(nn2(randomCells, query = sampledCell, k = 6))
    
    averageDistance = mean(neighbors$nn.dists[, 2:6])
    cv = CoefVar(as.matrix(neighbors$nn.dists[, 2:6]))
    
    averageDistances[j] = averageDistance
    cvs[j] = cv
    }
  
  averageDistances_subSims[ , i] = averageDistances
  cvs_subSims[ ,i] = cvs
  }

averageDistances_subSims = as.data.frame(averageDistances_subSims)
averageDistanceRandomSub = mean(unlist(averageDistances_subSims))
averageDistances_subSims = averageDistances_subSims %>% mutate_all(~ . -averageDistanceRandomSub) %>% mutate(typeA = "sim") %>% rename(value = V1)

cvs_subSims = as.data.frame(cvs_subSims)

avgDistances_singletsSub = avgDistances_singletsSub %>%  mutate(value = value- averageDistanceRandomSub)
avgDistances_bias = avgDistances_bias %>% mutate(value = value- averageDistanceRandomSub)

### combining all of the control data for plotting

allControls = rbind(avgDistances_singletsSub, avgDistances_bias, averageDistances_subSims)
allControls.summary = allControls %>% group_by(typeA) %>% summarise(sd = sd(as.numeric(value), na.rm = TRUE), value = mean(as.numeric(value), na.rm = TRUE)) #done above, this is maintained
allControls$typeA = factor(allControls$typeA, levels = c("bias", "singlets", "sim")) # ordering it in a more intuitive way for plotting

wilcox.test(avgDistances_bias$value, avgDistances_singletsSub$value, alternative = "less") #checking that biased NN distance is less than the singlet NN distance (p = 2.2e-16)
wilcox.test(avgDistances_singletsSub$value, averageDistances_subSims$value, alternative = "less") #checking if the singlet NN distance is less than the subsampled cell NN distance (p=0.9976)

plot = ggplot(allControls, aes(x = typeA, y = as.numeric(value))) +
  geom_point(aes(color = typeA), shape = 16, alpha = 0.7, size = 1, position = position_jitterdodge(jitter.width = 1.25, dodge.width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(color = "black", alpha = 0, size = 0.8) +
  theme_classic((base_size = 18)) +
  labs(x = "", y = "mean distance to 5NN") +
  geom_pointrange(aes(ymin = value - sd, ymax = value + sd, color = "black"), size = 1, color = "black", position = position_dodge(width = 0.5), data = allControls.summary) + #keeping data.summary for non-subsampled points!
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM01_sample3_biasControl.svg'), width = 4, height = 4) #### RENAME


### note average neighbor distance values can differ slightly based on the subsampling, without altering any of the conclusions.
