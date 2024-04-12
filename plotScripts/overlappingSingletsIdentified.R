### identifying which singlets are identified by which method for Zhang, Melzer et al 2024
### Created by Madeline E Melzer on 20240214
### Last edited by Madeline E Melzer on 20240319

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(ggtext) 

doubletListDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doublet_objects"
doubletInfoDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo2"
similarityResultsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo/similarity/"


##########################################################################################################################
# getting singlets for each method
##########################################################################################################################

method = "doublet_finder"
doubletFinderPath = file.path(doubletInfoDirectory, "doublet_finder")
files <- list.files(doubletFinderPath, pattern = "exp_0.1__act_0.1_", full.names = TRUE)

similarityResults <- data.frame(
  dataset = character(),
  sample = character(),
  lengthToCompare = integer(),
  labelDoublets = integer(),
  doubletFinderDoublets = integer(),
  scDblFinderDoublets = integer(),
  hybridDoublets = integer(),
  cxdsDoublets = integer(),
  bcdsDoublets = integer(),
  scrubletDoublets = integer(),
  minDoubletFinderScore = numeric(),
  minscDblFinderScore = numeric(),
  minHybridScore = numeric(),
  minCxdsScore = numeric(),
  minBcdsScore = numeric(),
  minScrubletScore = numeric(),
  similarityScore_doubletFinder_scDblFinder = numeric(),
  similarityScore_doubletFinder_hybrid = numeric(),
  similarityScore_doubletFinder_cxds = numeric(),
  similarityScore_doubletFinder_bcds = numeric(),
  similarityScore_doubletFinder_scrublet = numeric(),
  similarityScore_scDblFinder_hybrid = numeric(),
  similarityScore_scDblFinder_cxds = numeric(),
  similarityScore_scDblFinder_bcds = numeric(),
  similarityScore_scDblFinder_scrublet = numeric(),
  similarityScore_hybrid_cxds = numeric(),
  similarityScore_hybrid_bcds = numeric(),
  similarityScore_hybrid_scrublet = numeric(),
  similarityScore_cxds_bcds = numeric(),
  similarityScore_cxds_scrublet = numeric(),
  similarityScore_bcds_scrublet = numeric(),
  stringsAsFactors = FALSE
)


file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo2/doublet_finder/Biorxiv___1_DMSO_A___exp_0.1__act_0.1___doubletData.csv"
#file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo2/doublet_finder/watermelon___T47D-lag-1___exp_0.1__act_0.1___doubletData.csv"
#file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo2/doublet_finder/ClonMapper___FM1___exp_0.1__act_0.1___doubletData.csv"
#file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo2/doublet_finder/non_cancer___10kbarsiblingA___exp_0.1__act_0.1___doubletData.csv"  


for (file in files) {
  filename <- basename(file)
  
  # Extract dataset and sample names from the filename
  parts <- str_split(filename, "___", simplify = TRUE)
  dataset <- parts[1]
  sample <- parts[2]
  
  # DoubletFinder
  doubletFinderPath <- file.path("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo2/doublet_finder/")
  doublet_finder <- read.csv(paste0(doubletFinderPath, filename), stringsAsFactors = FALSE) %>%
    rename_with(~ if_else(str_detect(.x, "pANN"), "doubletFinderScore", .x), .cols = contains("pANN")) %>%
    rename_with(~ if_else(str_detect(.x, "DF.classifications"), "doubletFinderCall", .x), .cols = contains("DF.classifications")) %>%
    mutate(doubletFinderCall = tolower(doubletFinderCall))
  
  
  # scDblFinder
  scDblFinderPath <- file.path(doubletInfoDirectory, "scDblFinder/")
  scDblFinder <- read.csv(paste0(scDblFinderPath, filename), stringsAsFactors = FALSE) %>%
    rename_with(~ if_else(str_detect(.x, "scDblFinder.score"), "scDblFinderScore", .x), .cols = contains("scDblFinder.score")) %>%
    rename_with(~ if_else(str_detect(.x, "scDblFinder.class"), "scDblFinderCall", .x), .cols = contains("scDblFinder.class")) %>%
    mutate(scDblFinderCall = tolower(scDblFinderCall))
  
  
  # Hybrid, cxds, bcds
  hybridPath <- file.path(doubletInfoDirectory, "hybrid/")
  hybrid <- read.csv(paste0(hybridPath, filename), stringsAsFactors = FALSE) %>%
    rename_with(~ if_else(str_detect(.x, "hybrid_score"), "hybridScore", .x), .cols = contains("hybrid_score")) %>%
    rename_with(~ if_else(str_detect(.x, "hybrid_call"), "hybridCall", .x), .cols = contains("hybrid_call")) %>%
    rename_with(~ if_else(str_detect(.x, "cxds_score"), "cxdsScore", .x), .cols = contains("cxds_score")) %>%
    rename_with(~ if_else(str_detect(.x, "cxds_call"), "cxdsCall", .x), .cols = contains("cxds_call")) %>%
    rename_with(~ if_else(str_detect(.x, "bcds_score"), "bcdsScore", .x), .cols = contains("bcds_score")) %>%
    rename_with(~ if_else(str_detect(.x, "bcds_call"), "bcdsCall", .x), .cols = contains("bcds_call")) %>%
    mutate(hybridCall = tolower(hybridCall)) %>%
    mutate(bcdsCall = tolower(bcdsCall)) %>%
    mutate(cxdsCall = tolower(cxdsCall)) %>%
    mutate(across(c(hybridCall, bcdsCall, cxdsCall), ~ case_when(
      .x == "false" ~ "singlet",
      .x == "true" ~ "doublet",
      TRUE ~ .x  # Retain the original value if it's neither "true" nor "false"
    )))
  
  # Scrublet
  scrubletFilename <- str_replace(filename, "___doubletData", "")
  scrubletPath <- file.path("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo2/scrublet/act__0.1/")
  scrublet <- read.csv(paste0(scrubletPath, scrubletFilename), stringsAsFactors = FALSE) %>%
    rename_with(~ if_else(str_detect(.x, "score"), "scrubletScore", .x), .cols = contains("score")) %>%
    rename_with(~ if_else(str_detect(.x, "label"), "scrubletCall", .x), .cols = contains("label")) %>%
    rename_with(~ if_else(str_detect(.x, "barcode"), "cellID", .x), .cols = contains("barcode")) %>%
    mutate(scrubletCall = tolower(scrubletCall)) %>%
    mutate(across(c(scrubletCall), ~ case_when(
      .x == "false" ~ "singlet",
      .x == "true" ~ "doublet",
      TRUE ~ .x  # Retain the original value if it's neither "true" nor "false"
    )))
  
  first_join <- inner_join(doublet_finder, scDblFinder, by = c("cellID", "label", "nCount_RNA", "nFeature_RNA"))
  second_join <- inner_join(first_join, hybrid, by = c("cellID", "label", "nCount_RNA", "nFeature_RNA"))
  toCompare <- inner_join(second_join, scrublet, by = c("cellID"))
  
  # Calculate doublet similarity scores
  toCompare_doublets <- toCompare[toCompare$label == "doublet", ]
  
  doubletFinder_scDblFinder <- getSimilarityScore(toCompare_doublets, "doubletFinderCall", "scDblFinderCall")
  doubletFinder_hybrid <- getSimilarityScore(toCompare_doublets, "doubletFinderCall", "hybridCall")
  doubletFinder_cxds <- getSimilarityScore(toCompare_doublets, "doubletFinderCall", "cxdsCall")
  doubletFinder_bcds <- getSimilarityScore(toCompare_doublets, "doubletFinderCall", "bcdsCall")
  doubletFinder_scrublet <- getSimilarityScore(toCompare_doublets, "doubletFinderCall", "scrubletCall")
  scDblFinder_hybrid <- getSimilarityScore(toCompare_doublets, "scDblFinderCall", "hybridCall")
  scDblFinder_cxds <- getSimilarityScore(toCompare_doublets, "scDblFinderCall", "cxdsCall")
  scDblFinder_bcds <- getSimilarityScore(toCompare_doublets, "scDblFinderCall", "bcdsCall")
  scDblFinder_scrublet <- getSimilarityScore(toCompare_doublets, "scDblFinderCall", "scrubletCall")
  hybrid_cxds <- getSimilarityScore(toCompare_doublets, "hybridCall", "cxdsCall")
  hybrid_bcds <- getSimilarityScore(toCompare_doublets, "hybridCall", "bcdsCall")
  hybrid_scrublet <- getSimilarityScore(toCompare_doublets, "hybridCall", "scrubletCall")
  cxds_bcds <- getSimilarityScore(toCompare_doublets, "cxdsCall", "bcdsCall")
  cxds_scrublet <- getSimilarityScore(toCompare_doublets, "cxdsCall", "scrubletCall")
  bcds_scrublet <- getSimilarityScore(toCompare_doublets, "bcdsCall", "scrubletCall")
  
  # Calculate minimum scores for doublet calls
  minDoubletFinderScore <- min(doublet_finder$doubletFinderScore[doublet_finder$doubletFinderCall == "doublet"], na.rm = TRUE)
  minscDblFinderScore <- min(scDblFinder$scDblFinderScore[scDblFinder$scDblFinderCall == "doublet"], na.rm = TRUE)
  minHybridScore <- min(hybrid$hybridScore[hybrid$hybridCall == "doublet"], na.rm = TRUE)
  minCxdsScore <- min(hybrid$cxdsScore[hybrid$cxdsCall == "doublet"], na.rm = TRUE)
  minBcdsScore <- min(hybrid$bcdsScore[hybrid$bcdsCall == "doublet"], na.rm = TRUE)
  minScrubletScore <- min(scrublet$scrubletScore[scrublet$scrubletCall == "doublet"], na.rm = TRUE)
  
  # Count doublets
  labelDoublets <- sum(toCompare$label == "doublet", na.rm = TRUE)
  doubletFinderDoublets <- sum(toCompare$doubletFinderCall == "doublet", na.rm = TRUE)
  scDblFinderDoublets <- sum(toCompare$scDblFinderCall == "doublet", na.rm = TRUE)
  hybridDoublets <- sum(hybrid$hybridCall == "doublet", na.rm = TRUE)
  cxdsDoublets <- sum(hybrid$cxdsCall == "doublet", na.rm = TRUE)
  bcdsDoublets <- sum(hybrid$bcdsCall == "doublet", na.rm = TRUE)
  scrubletDoublets <- sum(scrublet$scrubletCall == "doublet", na.rm = TRUE)
  
  if(any(scrublet$scrubletCall == "doublet")) {
    min_value <- min(scrublet$scrubletScore[scrublet$scrubletCall == "doublet"], na.rm = TRUE)
    if(is.infinite(min_value)) {
      cat("No valid 'scrubletScore' for 'doublet' in file:", filename, "\n")
    }
  } else {
    # Print statement for identifying the file with no 'doublet' labels
    cat("No 'doublet' labels found in file:", filename, "\n")
  }
  
  
  # Inside your loop, after calculating all metrics for each file
  similarityResults <- rbind(similarityResults, data.frame(
    dataset = dataset,
    sample = sample,
    lengthToCompare = nrow(toCompare),
    labelDoublets = sum(toCompare$label == "doublet", na.rm = TRUE),
    doubletFinderDoublets = doubletFinderDoublets,
    scDblFinderDoublets = scDblFinderDoublets,
    hybridDoublets = hybridDoublets,
    cxdsDoublets = cxdsDoublets,
    bcdsDoublets = bcdsDoublets,
    scrubletDoublets = scrubletDoublets,
    minDoubletFinderScore = minDoubletFinderScore,
    minscDblFinderScore = minscDblFinderScore,
    minHybridScore = minHybridScore,
    minCxdsScore = minCxdsScore,
    minBcdsScore = minBcdsScore,
    minScrubletScore = minScrubletScore,
    similarityScore_doubletFinder_scDblFinder = doubletFinder_scDblFinder,
    similarityScore_doubletFinder_hybrid = doubletFinder_hybrid,
    similarityScore_doubletFinder_cxds = doubletFinder_cxds,
    similarityScore_doubletFinder_bcds = doubletFinder_bcds,
    similarityScore_doubletFinder_scrublet = doubletFinder_scrublet,
    similarityScore_scDblFinder_hybrid = scDblFinder_hybrid,
    similarityScore_scDblFinder_cxds = scDblFinder_cxds,
    similarityScore_scDblFinder_bcds = scDblFinder_bcds,
    similarityScore_scDblFinder_scrublet = scDblFinder_scrublet,
    similarityScore_hybrid_cxds = hybrid_cxds,
    similarityScore_hybrid_bcds = hybrid_bcds,
    similarityScore_hybrid_scrublet = hybrid_scrublet,
    similarityScore_cxds_bcds = cxds_bcds,
    similarityScore_cxds_scrublet = cxds_scrublet,
    similarityScore_bcds_scrublet = bcds_scrublet,
    stringsAsFactors = FALSE
  ))
  
  #print(similarityResults)
}


write_csv(similarityResults, paste0(similarityResultsDirectory, "similarityScoresAcrossDetectionMethods.csv"))


getSimilarityScore <- function(df, colName1, colName2) {
  matches <- df[[colName1]] == df[[colName2]]
  similarityScore <- sum(matches, na.rm = TRUE) / nrow(df)  # Ensure NA values are ignored
  return(similarityScore)
}


similarityResults_filtered = similarityResults %>% filter(labelDoublets !=0)


##########################################################################################################################
# plotting
##########################################################################################################################

results = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo/unfiltered/similaritysimilarityScoresAcrossDetectionMethods.csv")

results = results %>% #updated 20240321
  filter(
  !(
    (dataset == "SPLINTR" & sample %in% c("inVitro_KRAS", "retransplant")) |
      (dataset == "TREX" & sample == "brain1") |
      (dataset == "LARRY" & sample %in% c("d4_nBC", "d4_R_4")) |
      (dataset == "smartseq3" & sample %in% c("sample1"))
  ))

write.csv(results, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo/similarityScoresAcrossDetectionMethods_allSamples.csv")


pattern <- "dataset|similarityScore"

# Use grep to find columns that match the pattern
matching_columns <- grep(pattern, names(results), value = TRUE)

# Subset the dataframe to keep only the matching columns
df_filtered <- results[, matching_columns]

df_filtered <- df_filtered %>%
  group_by(dataset) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = 'drop')

write_csv(df_filtered, paste0(similarityResultsDirectory, "similarityScoresAcrossDetectionMethods_groupedByDataset.csv"))



methods <- unique(unlist(strsplit(sub("similarityScore_", "", matching_columns[-1]), "_")))



create_similarity_matrix <- function(dataset_name, results_df) {
  # Subset for the current dataset
  dataset_df <- results_df %>% filter(dataset == dataset_name)
  
  # Extract methods from column names
  methods <- colnames(dataset_df)
  methods <- methods[sapply(methods, function(x) grepl("similarityScore", x))]
  methods <- unique(unlist(strsplit(methods, "_")))
  methods <- methods[!methods %in% c("similarityScore", "dataset")]
  
  # Create an empty matrix
  similarity_matrix <- matrix(NA, length(methods), length(methods), 
                              dimnames = list(methods, methods))
  
  # Populate the matrix, taking the order into account
  for (i in 1:length(methods)) {
    for (j in i:length(methods)) {  # Start from i to avoid checking the order
      if (i != j) {
        # Create both possible column names, because we don't know the order
        score_col1 <- paste("similarityScore", methods[i], methods[j], sep = "_")
        score_col2 <- paste("similarityScore", methods[j], methods[i], sep = "_")
        
        # Check which column name exists, and retrieve the score
        if (score_col1 %in% colnames(dataset_df)) {
          similarity_matrix[i, j] <- dataset_df[[score_col1]]
        } else if (score_col2 %in% colnames(dataset_df)) {
          similarity_matrix[i, j] <- dataset_df[[score_col2]]
        } else {
          warning(paste("Neither", score_col1, "nor", score_col2, "exist in dataset", dataset_name))
        }
        
        # Assuming symmetry
        similarity_matrix[j, i] <- similarity_matrix[i, j]
      }
    }
  }
  
  # Set diagonal to 1, assuming perfect similarity with self
  diag(similarity_matrix) <- 1
  
  # Write the matrix to a CSV file
  write.csv(similarity_matrix, 
            file = paste0(similarityResultsDirectory, dataset_name, "_similarityMatrix.csv"), 
            row.names = TRUE)
  
  dataset_name
}

# Apply the function to each dataset
lapply(unique(df_filtered$dataset), create_similarity_matrix, results_df = df_filtered) #will return "NULL" because it doesnt return anything explicitly

##########################################################################################################################
plotting
##########################################################################################################################
plotsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/similarityScore/"

matrix <- read.csv(paste0(similarityResultsDirectory, "mean_similarityMatrix.csv"), row.names = 1)
dataset = "mean"

matrix <- as.matrix(matrix)

long_matrix = melt(matrix)

levels(long_matrix$Var1)[levels(long_matrix$Var1) == "scDblFinder"] <- "scDblFinder"
levels(long_matrix$Var2)[levels(long_matrix$Var2) == "scDblFinder"] <- "scDblFinder"
levels(long_matrix$Var1)[levels(long_matrix$Var1) == "doubletFinder"] <- "DoubletFinder"
levels(long_matrix$Var2)[levels(long_matrix$Var2) == "doubletFinder"] <- "DoubletFinder"
levels(long_matrix$Var1)[levels(long_matrix$Var1) == "scrublet"] <- "Scrublet"
levels(long_matrix$Var2)[levels(long_matrix$Var2) == "scrublet"] <- "Scrublet"

lower_triangle_data <- long_matrix %>%
  filter(as.numeric(Var2) >= as.numeric(Var1))

plot = ggplot(lower_triangle_data, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  geom_text(aes(label=sprintf("%.2f", value)), size=4, color="white", vjust=0.5, hjust=0.5) +
  scale_fill_gradient(low="#e4ffb5", high="#1a4600", limit = c(0.40, 1)) +
  theme_classic() +
  theme(axis.text.y = element_text(hjust=1)) +  # Rotate x axis labels if needed
  scale_x_discrete(position = "top") +
  labs(x='', y='', fill='Score')  # Customize axis labels and legend title
plot

ggsave(plot, file = paste0(plotsDirectory, dataset, "_similarityScore.svg"), width = 5, height = 4)
ggsave(plot, file = paste0(plotsDirectory, dataset, "_similarityScore.png"), width = 5, height = 4)











#### looping

# List all CSV files in the directory
csv_files <- list.files(similarityResultsDirectory, pattern = "\\.csv$", full.names = TRUE)

# Loop through the CSV files
for (file_path in csv_files) {
  # Extract the dataset name from the file name
  file_name <- basename(file_path)
  dataset <- sub("_similarityMatrix\\.csv$", "", file_name)
  
  # Read and process the CSV file
  matrix <- read.csv(file_path, row.names = 1)
  matrix <- as.matrix(matrix)
  long_matrix = melt(matrix)
  
  levels(long_matrix$Var1)[levels(long_matrix$Var1) == "scDblFinder"] <- "scDblFinder"
  levels(long_matrix$Var2)[levels(long_matrix$Var2) == "scDblFinder"] <- "scDblFinder"
  levels(long_matrix$Var1)[levels(long_matrix$Var1) == "doubletFinder"] <- "DoubletFinder"
  levels(long_matrix$Var2)[levels(long_matrix$Var2) == "doubletFinder"] <- "DoubletFinder"
  levels(long_matrix$Var1)[levels(long_matrix$Var1) == "scrublet"] <- "Scrublet"
  levels(long_matrix$Var2)[levels(long_matrix$Var2) == "scrublet"] <- "Scrublet"
  
  lower_triangle_data <- long_matrix %>%
    filter(as.numeric(Var2) >= as.numeric(Var1))
  
  # Create the heatmap plot
  plot = ggplot(lower_triangle_data, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    geom_text(aes(label=sprintf("%.2f", value)), size=4, color="white", vjust=0.5, hjust=0.5) +
    scale_fill_gradient(low="#e4ffb5", high="#1a4600", limit = c(0.35, 1)) +
    theme_classic() +
    theme(axis.text.y = element_text(hjust=1)) +  # Rotate x axis labels if needed
    scale_x_discrete(position = "top") +
    labs(x='', y='', fill='Score') + # Customize axis labels and legend title
    theme(legend.position = "none",
          axis.text.x = element_blank(),  # Remove x axis numbers
          axis.text.y = element_blank())
  plot
  print(plot)
  
  # Save the plots
  ggsave(plot, file = paste0(plotsDirectory, dataset, "_similarityHeatMap.svg"), width = 2, height = 2)
  ggsave(plot, file = paste0(plotsDirectory, dataset, "_similarityHeatMap.png"), width = 2, height = 2)
}

############## all datasets for main started 20240318, updated 20240321

allDatasets = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo/similarityScoresAcrossDetectionMethods_groupedByDataset.csv")

allDatasets_numeric = allDatasets[-1] #getting rid of "dataset" column.

## mean
mean = colMeans(allDatasets_numeric, na.rm = TRUE)
mean = t(mean)
mean = as_data_frame(mean)
names(mean) <- names(allDatasets_numeric)

## median
median_values <- apply(allDatasets_numeric, MARGIN = 2, FUN = median, na.rm = TRUE)
median <- as.data.frame(t(median_values))
names(median) <- names(allDatasets_numeric)


## i ran this for each median and mean

dataset_name = "mean"
dataset_df = mean

# Extract methods from column names
methods <- colnames(dataset_df)
methods <- methods[sapply(methods, function(x) grepl("similarityScore", x))]
methods <- unique(unlist(strsplit(methods, "_")))
methods <- methods[!methods %in% c("similarityScore", "dataset")]

# Create an empty matrix
similarity_matrix <- matrix(NA, length(methods), length(methods), 
                            dimnames = list(methods, methods))

# Populate the matrix, taking the order into account
for (i in 1:length(methods)) {
  for (j in i:length(methods)) {  # Start from i to avoid checking the order
    if (i != j) {
      # Create both possible column names, because we don't know the order
      score_col1 <- paste("similarityScore", methods[i], methods[j], sep = "_")
      score_col2 <- paste("similarityScore", methods[j], methods[i], sep = "_")
      
      # Check which column name exists, and retrieve the score
      if (score_col1 %in% colnames(dataset_df)) {
        similarity_matrix[i, j] <- dataset_df[[score_col1]]
      } else if (score_col2 %in% colnames(dataset_df)) {
        similarity_matrix[i, j] <- dataset_df[[score_col2]]
      } else {
        warning(paste("Neither", score_col1, "nor", score_col2, "exist in dataset", dataset_name))
      }
      
      # Assuming symmetry
      similarity_matrix[j, i] <- similarity_matrix[i, j]
    }
  }
}

# Set diagonal to 1, assuming perfect similarity with self
diag(similarity_matrix) <- 1

# Write the matrix to a CSV file
write.csv(similarity_matrix, 
          file = paste0(similarityResultsDirectory, dataset_name, "_similarityMatrix.csv"), 
          row.names = TRUE)


## getting average values for mean and median
mean_global = rowMeans(mean[, sapply(mean, is.numeric)], na.rm = TRUE)



#### plotting similarity score for each dataset in a column graph

plotsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/similarityScore/"

matrix <- read.csv(paste0(similarityResultsDirectory, "mean_similarityMatrix.csv"), row.names = 1)
matrix <- as.matrix(matrix)
long_matrix = melt(matrix)

long_matrix = long_matrix %>% filter(value != 1) %>% rename("method" = "Var1")
long_matrix_df = as_data_frame(long_matrix) %>%
  mutate(method = as.character(method)) %>%
  mutate(method = ifelse(method != "cxds", "combined", method))


summarized_df <- long_matrix_df %>%
  mutate(method = as.character(method)) %>%
  mutate(method = ifelse(method != "cxds", "combined", method)) %>%
  group_by(method) %>%
  summarize(mean_value = mean(value, na.rm = TRUE),
            std_dev = sd(value, na.rm = TRUE))

# Now 'summarized_df' contains the mean and standard deviation for each method

matrix_plot <- left_join(long_matrix_df, summarized_df, by = "method")

# Create the plot with point ranges and points
plot <- ggplot(matrix_plot, aes(x = method, y = value, group = method)) +
  geom_pointrange(aes(ymin = mean_value - std_dev, ymax = mean_value + std_dev, y = mean_value), shape = 16) +  # Point ranges for mean ± std_dev
  #geom_point(aes(y = mean_value), size = 3) +  # Points for the mean values
  geom_point(aes(y = value, color = method), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0)) +  # Dodge individual values with some transparency
  theme_classic() +
  labs(x = "", y = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotsDirectory, "similarityScore_averageDotPlot_cxdsAndCombined.svg"), width = 2, height = 3)
ggsave(plot, file = paste0(plotsDirectory, "similarityScore_averageDotPlot_cxdsAndCombined.png"), width = 2, height = 3)


#### for all points

allPoints = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo/similarityScoresAcrossDetectionMethods_allSamples.csv")

pattern <- "dataset|similarityScore"

# Use grep to find columns that match the pattern
matching_columns <- grep(pattern, names(allPoints), value = TRUE)

# Subset the dataframe to keep only the matching columns
df_filtered <- allPoints[, matching_columns]

allPoints_filtered = df_filtered %>% select(-dataset)


allPoint_filtered_long = melt(allPoints_filtered, variable.name = "methods", value.name = "value")

allPoint_filtered_long <- allPoint_filtered_long %>% #1365 rows = 91 samples * 15 method-method comparisons
  separate(methods, into = c("discard", "methodA", "methodB"), sep = "_") %>%
  select(-discard)

write.csv(allPoint_filtered_long, file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/doubletInfo/allPoints_similarityScores.csv") # 20240321

